use std::str::FromStr;

use anyhow::{anyhow, Result, bail};

// Partial implementation the GENCODE GTF file specification
// https://www.gencodegenes.org/pages/data_format.html
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
pub struct GTFRecord {
    // pub chromosome_name: String,
    // pub annotation_source: String,
    pub feature_type: FeatureType,
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
    // pub phase: Option<u8>,
    pub attributes: String
}

impl FromStr for GTFRecord {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        if s.starts_with('#') {
            bail!("Cannot parse comments, please remove them beforehand");
        }

        let mut cols = s.split('\t');

        let err: Result<Self> = Err(anyhow!("Syntax error, not valid GENCODE GTF: \"{s}\""));

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
                Some(s) => s.to_owned(),
                None => return err,
            },
        })
    }
}

impl GTFRecord {
    pub fn len(&self) -> usize {
        (self.start as i64 - self.end as i64).abs().try_into().unwrap()
    }
}