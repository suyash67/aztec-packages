//! Baseline numbers for the production Zcash Orchard prover (orchard 0.16.0 / halo2_proofs 0.4).
//!
//! Builds bundles with exactly `n` actions (one output per action, dummy spends, no padding), then times
//! `create_proof` and single-proof `verify` with the post-NU6.2 (`FixedPostNu6_2`) circuit, exactly as
//! `orchard/benches/circuit.rs` does, and records the serialized proof length.
//!
//! Usage: orchard_bench [comma-separated action counts] > results.json

use std::time::Instant;

use orchard::{
    builder::{Builder, BundleType},
    bundle::BundleVersion,
    circuit::{OrchardCircuitVersion, ProvingKey, VerifyingKey},
    keys::{FullViewingKey, Scope, SpendingKey},
    value::NoteValue,
    Anchor, Bundle,
};
use rand::{rand_core::UnwrapErr, rngs::SysRng};

fn median(mut xs: Vec<f64>) -> f64 {
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    xs[xs.len() / 2]
}

fn main() {
    let counts: Vec<usize> = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "1,2,4,8,16".to_string())
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect();
    let threads = std::env::var("RAYON_NUM_THREADS").unwrap_or_else(|_| "default".to_string());
    let mut rng = UnwrapErr(SysRng);

    let t = Instant::now();
    let vk = VerifyingKey::build(OrchardCircuitVersion::FixedPostNu6_2);
    let pk = ProvingKey::build(OrchardCircuitVersion::FixedPostNu6_2);
    let keygen_s = t.elapsed().as_secs_f64();
    eprintln!("keygen: {keygen_s:.3}s");

    let sk = SpendingKey::from_bytes([7; 32]).unwrap();
    let recipient = FullViewingKey::from(&sk).address_at(0u32, Scope::External);

    let mut rows = Vec::new();
    for &n in &counts {
        let mut builder = Builder::new(
            BundleType::Transactional {
                bundle_required: true,
                pad_to_minimum: Some(1),
            },
            BundleVersion::orchard_v2(),
            BundleVersion::orchard_v2().default_flags(),
            Anchor::from_bytes([0; 32]).unwrap(),
        )
        .unwrap();
        for _ in 0..n {
            builder
                .add_output(None, recipient, NoteValue::from_raw(10), [0; 512])
                .unwrap();
        }
        let bundle: Bundle<_, i64> = builder.build(&mut rng).unwrap().unwrap().0;
        assert_eq!(bundle.actions().len(), n);
        let instances: Vec<_> = bundle
            .actions()
            .iter()
            .map(|a| a.to_instance(*bundle.flags(), *bundle.anchor()))
            .collect();

        let reps = if n <= 4 { 5 } else { 3 };
        let mut prove_times = Vec::new();
        let mut proof = None;
        for _ in 0..reps {
            let t = Instant::now();
            let p = bundle
                .authorization()
                .create_proof(&pk, &instances, &mut rng)
                .unwrap();
            prove_times.push(t.elapsed().as_secs_f64());
            proof = Some(p);
        }
        let proof = proof.unwrap();
        let proof_bytes = proof.as_ref().len();

        let mut verify_times = Vec::new();
        for _ in 0..20 {
            let t = Instant::now();
            proof.verify(&vk, &instances).expect("valid proof");
            verify_times.push(t.elapsed().as_secs_f64());
        }

        let row = serde_json::json!({
            "system": "zcash-orchard-halo2",
            "actions": n,
            "threads": threads,
            "prove_s": median(prove_times.clone()),
            "prove_all_s": prove_times,
            "verify_s": median(verify_times.clone()),
            "proof_bytes": proof_bytes,
            "keygen_s": keygen_s,
        });
        eprintln!("{row}");
        rows.push(row);
    }
    println!("{}", serde_json::to_string_pretty(&rows).unwrap());
}
