// Partial implementation the GENCODE GTF file specification
// https://www.gencodegenes.org/pages/data_format.html

use anyhow::{anyhow, bail, Result};
use flate2::read::GzDecoder;
use lazy_static::lazy_static;
use regex::Regex;
use std::str::FromStr;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

#[derive(Debug, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

impl FromStr for Strand {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Plus),
            "-" => Ok(Strand::Minus),
            _ => todo!(),
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum FeatureType {
    Gene,
    Transcript,
    Exon,
    CDS,
    UTR,
    StartCodon,
    StopCodon,
    Selenocysteine,
}

impl FromStr for FeatureType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use FeatureType::*;

        match s {
            "gene" => Ok(Gene),
            "transcript" => Ok(Transcript),
            "exon" => Ok(Exon),
            "CDS" => Ok(CDS),
            "UTR" => Ok(UTR),
            "start_codon" => Ok(StartCodon),
            "stop_codon" => Ok(StopCodon),
            "Selenocysteine" => Ok(Selenocysteine),
            _ => Err(anyhow!("Couldn't parse feature type \"{s}\"")),
        }
    }
}

#[derive(Debug)]
pub struct Attributes {
    pub gene_id: String,
    pub transcript_id: String,
}

impl FromStr for Attributes {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        lazy_static! {
            static ref RE: Regex =
                Regex::new("gene_id \"(ENSG[^\"]+)\";.*transcript_id \"(ENST[^\"]+)\";").unwrap();
        }

        let captures = RE.captures(s).unwrap();

        Ok(Self {
            gene_id: captures[1].to_string(),
            transcript_id: captures[2].to_string(),
        })
    }
}

#[derive(Debug)]
pub struct GTFRecord {
    // pub chromosome_name: String,
    // pub annotation_source: String,
    pub feature_type: FeatureType,
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
    // pub phase: Option<u8>,
    pub attributes: Attributes,
}

impl FromStr for GTFRecord {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        if s.starts_with('#') {
            bail!("Cannot parse comments, please remove them beforehand");
        }

        let mut cols = s.split('\t');

        let err: Result<Self> = Err(anyhow!("Syntax error, not valid GENCODE GTF: \"{s}\""));

        Ok(Self {
            feature_type: match cols.nth(2) {
                Some(s) => s.parse::<FeatureType>()?,
                None => return err,
            },
            start: match cols.next() {
                Some(s) => s.parse::<usize>()?,
                None => return err,
            },
            end: match cols.next() {
                Some(s) => s.parse::<usize>()?,
                None => return err,
            },
            strand: match cols.nth(1) {
                Some(s) => s.parse::<Strand>()?,
                None => return err,
            },
            attributes: match cols.nth(1) {
                Some(s) => s.parse::<Attributes>()?,
                None => return err,
            },
        })
    }
}

impl GTFRecord {
    pub fn len(&self) -> usize {
        (self.start as i64 - self.end as i64)
            .abs()
            .try_into()
            .unwrap()
    }
}

pub fn read_gtf_file(annotations: &Path) -> Result<Vec<GTFRecord>> {
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

pub fn get_longest_transcripts(gtf_records: &[GTFRecord]) -> Result<Vec<&GTFRecord>> {
    let mut longest_transcript_candidates: HashMap<&str, &GTFRecord> = HashMap::new(); // (gene_id, longest_transcript_so_far)

    let transcripts = gtf_records
        .iter()
        .filter(|&r| r.feature_type == FeatureType::Transcript);

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
