mod gc_content;
mod gtf;
mod plots;

use anyhow::{anyhow, Result};
use bio::io::fasta;
use clap::Parser;
use flate2::read::GzDecoder;
use regex::Regex;
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use crate::gtf::GTFRecord;

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

    let fasta_record = get_fasta_record(&cli.sequence)?;

    if cli.total_gc {
        plots::plot_gc_content(&fasta_record, 5000)?;
    }

    if cli.promotor_gc {
        get_promotor_regions(&cli.annotations)?;
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

fn get_promotor_regions(annotations: &Path) -> anyhow::Result<()> {
    let reader = MaybeCompressedReader::new(annotations)?;

    let regex = Regex::new(r"^chr1\t\w*\t(gene|transcript|start_codon)").unwrap();

    let records: Result<Vec<GTFRecord>> = reader
        .lines()
        .map(|line| line.unwrap())
        .filter(|line| regex.is_match(line))
        .map(|line| line.parse::<GTFRecord>())
        .collect();
    
    let records = records?;

    println!("{:#?}", records[0]);

    Ok(())
}

/// Reader for a file that could be gzip compressed or not
struct MaybeCompressedReader;

impl MaybeCompressedReader {
    /// Checks if the file is gz compressed by looking at the extension,
    /// and returns the correct Reader
    fn new(path: &Path) -> Result<Box<dyn BufRead>> {
        println!("Extension {:?}", path.extension().unwrap());
        match path.extension() {
            Some(ext) if ext == "gz" => {
                Ok(Box::new(BufReader::new(GzDecoder::new(File::open(path)?))))
            }
            _ => Ok(Box::new(BufReader::new(File::open(path)?))),
        }
    }
}