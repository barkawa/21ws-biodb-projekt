mod gc_content;
mod gtf;
mod plots;

use anyhow::{anyhow, Result};
use bio::io::fasta;
use clap::Parser;
use flate2::read::GzDecoder;
use regex::Regex;
use std::{
    collections::{HashMap, HashSet},
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

    let fasta_record = read_fasta_record(&cli.sequence)?;

    if cli.total_gc {
        plots::plot_gc_content(&fasta_record, 5000)?;
    }

    if cli.promotor_gc {
        let records = read_gtf_file(&cli.annotations)?;
        let longest_transcripts = get_longest_transcripts(&records)?;

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
    }

    if cli.promotor_nsome_affinity {
        todo!();
    }

    if cli.tfbs_nsome_affinity {
        todo!();
    }

    Ok(())
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

fn read_gtf_file(annotations: &Path) -> Result<Vec<GTFRecord>> {
    let reader = MaybeCompressedReader::new(annotations)?;

    // Select the lines we need with a regex first, and parse later (for performance reasons)
    let regex = Regex::new(
        r"^chr1\t(?:HAVANA|ENSEMBL)\t(?:transcript|start_codon).*gene_type..protein_coding.;",
    )
    .unwrap();

    let records: Result<Vec<GTFRecord>> = reader
        .lines()
        .map(|line| line.unwrap())
        .filter(|line| !line.starts_with('#'))
        .filter(|line| regex.is_match(line))
        .map(|line| line.parse::<GTFRecord>())
        .collect();

    records
}

fn get_longest_transcripts(gtf_records: &[GTFRecord]) -> Result<Vec<&GTFRecord>> {
    let mut longest_transcript_candidates: HashMap<&str, &GTFRecord> = HashMap::new();

    let transcripts = gtf_records
        .iter()
        .filter(|&r| r.feature_type == gtf::FeatureType::Transcript);

    for candidate in transcripts {
        longest_transcript_candidates
            .entry(&candidate.attributes.gene_id)
            .and_modify(|current| {
                if candidate.len() > current.len() {
                    *current = candidate;
                }
            })
            .or_insert(candidate);
    }

    let longest_transcripts: Vec<_> = longest_transcript_candidates.values().map(|x| *x).collect();

    Ok(longest_transcripts)
}

/// Reader for a file that could be gzip compressed or not
struct MaybeCompressedReader;

impl MaybeCompressedReader {
    /// Checks if the file is gz compressed by looking at the extension,
    /// and returns the correct Reader
    fn new(path: &Path) -> Result<Box<dyn BufRead>> {
        match path.extension() {
            Some(ext) if ext == "gz" => {
                Ok(Box::new(BufReader::new(GzDecoder::new(File::open(path)?))))
            }
            _ => Ok(Box::new(BufReader::new(File::open(path)?))),
        }
    }
}
