use std::str::FromStr;

use anyhow::{anyhow, Result};

// Partial implementation the GENCODE GTF file specification
// https://www.gencodegenes.org/pages/data_format.html
#[derive(Debug)]
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

#[derive(Debug)]
pub enum FeatureType {
    Gene,
    Transcript,
    // Exon,
    // CDS,
    // UTR,
    StartCodon,
    // StopCodon,
    // Selenocysteine,
}

impl FromStr for FeatureType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "gene" => Ok(FeatureType::Gene),
            "transcript" => Ok(FeatureType::Transcript),
            "start_codon" => Ok(FeatureType::StartCodon),
            _ => todo!(),
        }
    }
}

#[derive(Debug)]
pub struct GTFRecord {
    // chromosome_name: String,
    // annotation_source: String,
    feature_type: FeatureType,
    start: usize,
    end: usize,
    strand: Strand,
    // phase: Option<u8>,
    attributes: String
}

impl FromStr for GTFRecord {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let mut cols = s.split('\t');

        let err: Result<Self> = Err(anyhow!("Error parsing GTF record"));

        Ok(GTFRecord {
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
                //Some(s) => s.split(';').map(|s| s.to_string()).collect(),
                Some(s) => s.to_owned(),
                None => return err,
            },
        })
    }
}