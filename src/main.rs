use rand::Rng;
use rand::distr::{Distribution, Uniform};
use rand::prelude::IndexedRandom; 
use std::ops::Range;
use std::io::Write; 
use std::fs::File; 
use std::error::Error;

use clap::Parser;

struct PairedFastqRecord {
    id: String,
    seq: Vec<u8>,
    mate: Vec<u8>,
    quality_1: Vec<u8>,
    quality_2: Vec<u8>,
}

#[derive(Debug)]
struct FastqGenerator {
    bases: &'static [u8],
    read_length: usize,
    quality_range: Range<u8>
}

impl FastqGenerator {
    fn new(read_length: usize) -> Self {
        let phred_range: Range<u8> = 33u8..74u8;
        FastqGenerator { 
            bases: b"ATCG", 
            read_length, 
            quality_range: phred_range 
        }
    }

    fn sample_quality(&self, rng: &mut impl Rng) -> Vec<u8> {
        let dist = Uniform::new(self.quality_range.start, self.quality_range.end).unwrap();

        (0..self.read_length)
            .map(|_| dist.sample(rng)) 
            .collect()
    }

    fn sample_seq(&self, rng: &mut impl Rng) -> Vec<u8> {
        (0..self.read_length)
        .map(|_| {
            // Rand 0.9: choose() now requires the IndexedRandom trait import
            *self.bases.choose(rng).unwrap()
        })
        .collect()
    }

    fn generate_paired_record(&self, rng: &mut impl Rng, id_index: i32) -> PairedFastqRecord {
        let seq = self.sample_seq(rng);
        let qual_1 = self.sample_quality(rng);

        let mate = reverse_complement(&seq);
        
        // The quality for R2 is reversed from R1
        let qual_2: Vec<u8> = qual_1.iter().rev().copied().collect();

        PairedFastqRecord { 
            id: format!("READ_{:06}", id_index), 
            seq, 
            mate, 
            quality_1: qual_1,
            quality_2: qual_2
        }
    }
}

// --- Sequence Helper Functions ---

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => base
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev()
    .map(|base| complement(*base))
    .collect()
}

#[derive(Parser, Debug)]
#[command(version, about = "A simple tool to generate random paired-end fastq files.", long_about = None)]
struct Args {
    #[arg(short, long, default_value_t = String::from("synthetic_reads"))]
    outfile: String,
    
    #[arg(short, default_value_t = 1000)]
    n: i32,
    
    #[arg(short = 'l', default_value_t = 150)]
    read_len: i32
}

fn main() -> Result<(), Box<dyn Error>> { 
    let args = Args::parse();

    let output_file_prefix = args.outfile;
    let num_reads = args.n;
    let read_length = args.read_len;

    if num_reads <= 0 || read_length <= 0 {
        return Err("Number of reads and read length must be positive.".into());
    }

    let read_length_usize = read_length as usize;
    let num_reads_i32 = num_reads; 

    let generator = FastqGenerator::new(read_length_usize);

    // Rand 0.9: thread_rng() is renamed to rng()
    let mut rng = rand::rng();
    
    let r1_filepath = format!("{}_R1.fastq", output_file_prefix);
    let r2_filepath = format!("{}_R2.fastq", output_file_prefix);

    let mut r1_file = File::create(&r1_filepath)?;
    let mut r2_file = File::create(&r2_filepath)?;

    println!("Starting generation of {} paired reads (Length: {})", num_reads, read_length);

    for i in 0..num_reads_i32 {
        let record = generator.generate_paired_record(&mut rng, i);

        // --- Write Read 1 ---
        r1_file.write_all(format!("@{} /1\n", record.id).as_bytes())?;
        r1_file.write_all(&record.seq)?;
        r1_file.write_all(b"\n+\n")?; // Combined newlines for speed
        r1_file.write_all(&record.quality_1)?;
        r1_file.write_all(b"\n")?;

        // --- Write Read 2 ---
        r2_file.write_all(format!("@{} /2\n", record.id).as_bytes())?; 
        r2_file.write_all(&record.mate)?;
        r2_file.write_all(b"\n+\n")?;
        r2_file.write_all(&record.quality_2)?;
        r2_file.write_all(b"\n")?;
    }

    println!("🦀 Wrote {} paired reads of length {} to {}_R[12].fastq", num_reads, read_length, output_file_prefix);

    Ok(())
}