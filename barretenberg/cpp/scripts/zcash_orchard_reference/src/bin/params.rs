// Prints halo2's IPA parameters (k = 4) for Vesta: G_0, G_1, G_15, W, U as compressed points (hex, little-endian).
use halo2_proofs::poly::commitment::Params;
use pasta_curves::vesta;

fn main() {
    let params: Params<vesta::Affine> = Params::new(4);
    let mut buf = vec![];
    params.write(&mut buf).unwrap();
    // layout: k (u32 LE), then n generators, then n lagrange generators, then w, u (32 bytes each)
    let n = 16usize;
    let pt = |i: usize| hex(&buf[4 + 32 * i..4 + 32 * (i + 1)]);
    println!("G0 {}", pt(0));
    println!("G1 {}", pt(1));
    println!("G15 {}", pt(15));
    println!("W {}", pt(2 * n));
    println!("U {}", pt(2 * n + 1));
}

fn hex(b: &[u8]) -> String {
    b.iter().map(|x| format!("{:02x}", x)).collect()
}
