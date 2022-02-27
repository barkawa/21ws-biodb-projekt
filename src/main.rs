mod gc_content;
mod gtf;
mod plots;

use anyhow::{anyhow, Result};
use bio::io::fasta;
use clap::Parser;
use std::{
    collections::HashSet,
    fs::File,
    path::{Path, PathBuf},
};

#[derive(Parser)]
#[clap(author)]
struct Cli {
    /// File contatining a sequence for a single chromosome, in FASTA format (.fa/.fasta)
    sequence: PathBuf,

    /// File contatining sequence annotations in (gz compressed) GTF format (.gtf/.gtf.gz)
    annotations: PathBuf,

    /// File containing MNase-seq data in Wig format (.wig)
    mnase_seq: PathBuf,

    /// Plot the GC content for the whole chromosome
    #[clap(long)]
    total_gc: bool,

    /// Plot the average GC content of all promotor regions
    #[clap(long)]
    promotor_gc: bool,

    /// Plot the average Nucleosome affinity of of all promotor regions
    #[clap(long)]
    promotor_nsome_affinity: bool,

    /// Plot the average Nucleosome affinity of the TFBS
    #[clap(long)]
    tfbs_nsome_affinity: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let fasta_record = read_fasta_record(&cli.sequence)?;

    if cli.total_gc {
        plots::plot_gc_content(&fasta_record, 5000)?;
    }

    let records = gtf::read_gtf_file(&cli.annotations)?;
    let longest_transcripts = gtf::get_longest_transcripts(&records)?;

    let ids: HashSet<_> = HashSet::from_iter(
        longest_transcripts
            .iter()
            .map(|&r| &r.attributes.transcript_id),
    );

    let start_codons: Vec<_> = records
        .iter()
        .filter(|&r| r.feature_type == gtf::FeatureType::StartCodon)
        .filter(|&r| ids.contains(&r.attributes.transcript_id))
        .collect();

    let mut raw_promotors: Vec<PromotorRegion> = Vec::with_capacity(start_codons.len());

    // ATG is at 1000, CAT at 98
    // in fasta: ATG at p.start+1000, CAT at p.end+98
    for start_codon in &start_codons {
        if start_codon.strand == gtf::Strand::Plus {
            let start = start_codon.start - 1001;
            let end = start_codon.start + 100;

            raw_promotors.push(PromotorRegion {
                sequence: &fasta_record.seq()[start..end],
                location: start,
                strand: gtf::Strand::Plus,
            });
        } else {
            let start = start_codon.end + 1000;
            let end = start_codon.end - 101;

            raw_promotors.push(PromotorRegion {
                sequence: &fasta_record.seq()[end..start],
                location: end,
                strand: gtf::Strand::Minus,
            });
        }
    }

    if cli.promotor_gc {
        let promotors = raw_promotors
            .iter()
            .map(|p| match p.strand {
                Plus => p.sequence,
                Minus => &p.get_opposite_sequence(),
            })
            .map(|p| gc_content::get_gc_content(p, 150, 1));
            

        // for p in &promotors {
        //     println!(
        //         "{} {}",
        //         String::from_utf8_lossy(&p.sequence[1000..1003]),
        //         if p.start < p.end {
        //             String::from_utf8_lossy(&fasta_record.seq()[p.start+1000..p.end-98])
        //         } else {
        //             String::from_utf8_lossy(&fasta_record.seq()[p.end+98..p.start-1000])
        //         }
        //     );
        // }
    }

    if cli.promotor_nsome_affinity {
        todo!();
    }

    if cli.tfbs_nsome_affinity {
        todo!();
    }

    Ok(())
}

#[derive(Debug)]
struct PromotorRegion<'a> {
    sequence: &'a [u8],
    location: usize,
    strand: gtf::Strand,
}

impl PromotorRegion<'_> {
    fn get_opposite_sequence(&self) -> Vec<u8> {
        self.sequence
            .into_iter()
            .rev()
            .map(|&x| match x {
                b'A' => b'T',
                b'T' => b'A',
                b'G' => b'C',
                b'C' => b'G',
                other => other,
            })
            .collect()
    }
}

fn read_fasta_record(path: &Path) -> Result<fasta::Record> {
    let fasta_file = File::open(path)?;
    let mut fasta_records = fasta::Reader::new(fasta_file).records();

    match fasta_records.next() {
        Some(rec) => {
            if fasta_records.next().is_none() {
                Ok(rec?)
            } else {
                Err(anyhow!("Error: FASTA file has more than one record"))
            }
        }
        None => Err(anyhow!("Error: FASTA file has no records")),
    }
}
