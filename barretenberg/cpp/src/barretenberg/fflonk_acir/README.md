# fflonk_acir: Noir onto the three-wire arithmetization

`bb prove --scheme fflonk` needs an ACIR frontend for `../fflonk/`, which has three wires and five
selectors. Writing a second frontend would mean a second place for the frontend to be wrong, so this
module reuses the one that exists: ACIR builds an `UltraCircuitBuilder`, and `lower()` rewrites its
trace.

## What it does

Each Ultra arithmetic row spans five wire slots,

```
q_m w_l w_r + q_l w_l + q_r w_r + q_o w_o + q_4 w_4 + q_c [+ w_4(next row)] = 0
```

which three wires cannot hold, so it becomes two three-wire rows (or three, when the gate reaches
into the next row) threaded by fresh intermediates. Copy constraints carry over by construction: the
map from Ultra variables to three-wire variables is keyed on `real_variable_index`, the representative
after Ultra's own equality merging, so two slots Ultra identifies land on one variable here.

The lowered circuit is run through `Trace::check` before it is proven, so a mistranslated selector or
a dropped copy constraint surfaces as a named circuit violation rather than as a proof that fails
with no explanation.

Arithmetic gates and Poseidon2 rounds are covered. A Poseidon2 round becomes an S-box chain per lane
(`x^5` is three multiplications) and then the matrix, about thirty three-wire rows per Ultra row.

## What it refuses, and why each is not a gap waiting to be filled

**Lookups** — bitwise operations, SHA-256, Blake2s, Keccak, Pedersen. There is no three-wire encoding
of a lookup that is not astronomically expensive: expanding one into explicit constraints costs the
size of the table.

**Range constraints** — and this one is subtler than it looks, because the delta-range *quartic*
expands into three wires without trouble. The quartic is only half of Ultra's range argument. The
other half is that the sorted list the deltas run over is a permutation of the constrained variables,
which Ultra enforces with `real_variable_tags` and `tau`: a multiset equality layered onto sigma,
entirely separate from the copy constraints in `real_variable_index`.

This pass carries `real_variable_index` and nothing else, so lowering only the quartic would leave a
prover free to put any sorted list it liked in those wires — a soundness break, silent, and
indistinguishable from working. It is refused instead. Carrying the tags would need the three-wire
builder to express a permutation cycle that is not an equality class, which its union-find cannot;
that is an arithmetization change to `../fflonk/`, not a gap in this file.

**RAM/ROM** — dynamically indexed arrays need the same kind of argument, and the ROM-LogUp part is a
trace-sum subrelation the three-wire system has no mechanism for at all.

**Elliptic and non-native field gates** — expandable in principle, not yet done.

**A gate reading the next trace row at a block boundary**, where which row comes next is a property of
the trace layout rather than of the block, and is refused rather than guessed at.

So this backend serves programs built from field arithmetic and Poseidon2 that want the smallest
possible on-chain verifier — 207,189 execution gas, measured — and the smallest possible trusted code
base. **Integer types and comparisons compile to range constraints, so most real Noir programs belong
on `../ultra_fflonk/`**, which costs about 75k more execution gas and handles everything.
