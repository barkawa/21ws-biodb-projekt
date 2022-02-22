mod gc_content;

use bio::io::fasta;
use clap::{Parser, Subcommand};
use std::{fs::File, path::PathBuf};

#[derive(Parser)]
#[clap(author)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
    // /// File contatining a sequence for a single chromosome, in FASTA format (.fa/.fasta)
    // #[clap(short, long)]
    // sequence: PathBuf,

    // /// File contatining sequence annotations in compressed GTF format (.gtf.gz)
    // #[clap(short, long)]
    // annotations: PathBuf,

    // /// File containing MNase-seq data in Wig format (.wig)
    // #[clap(short, long)]
    // mnase_seq: PathBuf,
}

#[derive(Subcommand)]
enum Commands {
    /// Plots the GC content of a sequence
    GCPlot { file: PathBuf }
}

fn main() {
    let cli = Cli::parse();

    let fasta_file = match &cli.command {
        Commands::GCPlot { file } => file,
    };

    let fasta_file = File::open(fasta_file).unwrap();
    let mut fasta_records = fasta::Reader::new(fasta_file).records();
    let chr1 = fasta_records.next().unwrap().unwrap();

    // we only want one fasta record here
    assert!(fasta_records.next().is_none());

    match gc_content::plot_gc_content(&chr1, 5000) {
        Ok(()) => (),
        Err(why) => eprintln!("Plotting GC-Content failed: {why}"),
    }
}
