mod gc_content;
mod plots;

use anyhow::{bail, Result};
use bio::io::fasta;
use clap::Parser;
use std::{fs::File, path::{PathBuf, Path}};

#[derive(Parser)]
#[clap(author)]
struct Cli {
    /// File contatining a sequence for a single chromosome, in FASTA format (.fa/.fasta)
    sequence: PathBuf,

    /// File contatining sequence annotations in gz compressed GTF format (.gtf.gz)
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

fn main() -> Result<()>{
    let cli = Cli::parse();

    let record = get_fasta_record(&cli.sequence)?;

    if cli.total_gc {
        plots::plot_gc_content(&record, 5000)?;
    }

    if cli.promotor_gc {
        todo!();
    }

    if cli.promotor_nsome_affinity {
        todo!();
    }

    if cli.tfbs_nsome_affinity {
        todo!();
    }

    Ok(())
}

fn get_fasta_record(path: &Path) -> Result<fasta::Record> {
    let fasta_file = File::open(path)?;
    let mut fasta_records = fasta::Reader::new(fasta_file).records();

    let record = match fasta_records.next() {
        Some(rec) => {
            if fasta_records.next().is_some() {
                bail!("Error: FASTA file has more than one record")
            }
            rec?
        },
        None => bail!("Error: FASTA file has no records"),
    };

    Ok(record)
}