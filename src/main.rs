mod gc_content;
mod gtf;
mod plots;

use anyhow::{anyhow, Result};
use bigtools::{bigwigread::BigWigRead, seekableread::ReopenableFile};
use bio::io::fasta;
use clap::Parser;
use regex::Regex;
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

    /// File containing MNase-seq data in BigWig format (.bigWig)
    mnase_seq: String,

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
        plots::plot1(&fasta_record, 5000)?;
        return Ok(());
    }

    let records = gtf::read_gtf_file(&cli.annotations)?;
    let longest_transcripts = gtf::get_longest_transcripts(&records)?;
    let promotors = find_promotor_regions(&fasta_record, &records, longest_transcripts);

    if cli.promotor_gc {
        let mut total_gc: Vec<(usize, f64)> = vec![(0, 0.0); promotors[0].sequence.len() - 149];

        for p in &promotors {
            let sequence = if p.strand == gtf::Strand::Minus {
                p.get_opposite_sequence()
            } else {
                p.sequence.clone()
            };

            let gc: Vec<_> = gc_content::get_gc_content(&sequence, 150, 1).collect();

            for ((i, x), (j, y)) in total_gc.iter_mut().zip(gc.iter()) {
                *i = *j;
                *x += y;
            }
        }

        let total_gc: Vec<_> = total_gc
            .into_iter()
            .map(|(i, x)| (i, x / promotors.len() as f64))
            .collect();

        plots::plot2(&total_gc)?;
        return Ok(());
    }

    let mut big_wig_reader = BigWigRead::from_file_and_attach(cli.mnase_seq.as_str()).unwrap();

    if cli.promotor_nsome_affinity {
        let mut total_affinity = [(0u32, 0.0); 1100];

        for p in &promotors {
            let mut affinity =
                big_wig_reader.values("chr1", p.location as u32, p.location as u32 + 1100)?;

            if p.strand == gtf::Strand::Minus {
                affinity.reverse();
            }

            for (i, v) in affinity.iter().enumerate() {
                if !v.is_nan() {
                    total_affinity[i].1 += *v as f64;
                    total_affinity[i].0 += 1;
                }
            }
        }

        let avg_affinity =
            (-1000..100).zip(total_affinity.iter().map(|&(count, v)| v / count as f64));

        plots::plot3(avg_affinity)?;
        return Ok(());
    }

    if cli.tfbs_nsome_affinity {
        let regex_ap_1 = regex::bytes::Regex::new(r"TGA(C|G)TCA").unwrap();
        let regex_nf_y = regex::bytes::Regex::new(r"(CCAAT|ATTGG)").unwrap();

        let ap_1 = get_tfbs_avg_nsome_affinity(&regex_ap_1, &promotors, &mut big_wig_reader)?;
        plots::plot4(ap_1, "AP-1")?;

        let nf_y = get_tfbs_avg_nsome_affinity(&regex_nf_y, &promotors, &mut big_wig_reader)?;
        plots::plot4(nf_y, "NF-Y")?;
    }

    Ok(())
}

fn get_tfbs_avg_nsome_affinity(
    regex: &regex::bytes::Regex,
    promotors: &[PromotorRegion],
    big_wig_reader: &mut BigWigRead<ReopenableFile, File>,
) -> Result<Vec<(i32, f64)>> {
    let mut total_affinity = vec![(0, 0.0); 1001];

    for p in promotors {
        for m in regex.find_iter(&p.sequence) {
            let tfbs_center = p.location + ((m.end() - m.start()) / 2);

            let mut affinity = big_wig_reader.values(
                "chr1",
                tfbs_center as u32 - 500,
                tfbs_center as u32 + 501,
            )?;

            // if p.strand == gtf::Strand::Minus {
            //     affinity.reverse();
            // }

            for (i, v) in affinity.iter().enumerate() {
                if !v.is_nan() {
                    total_affinity[i].1 += *v as f64;
                    total_affinity[i].0 += 1;
                }
            }
        }
    }
    
    println!("{total_affinity:#?}");

    let avg_affinity = (-500..501).zip(total_affinity.iter().map(|&(count, v)| v / count as f64)).collect();

    Ok(avg_affinity)
}

#[derive(Debug)]
struct PromotorRegion {
    sequence: Vec<u8>,
    location: usize, // first index in fasta file, ignoring direction
    strand: gtf::Strand,
}

impl PromotorRegion {
    fn get_opposite_sequence(&self) -> Vec<u8> {
        self.sequence
            .iter()
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

fn find_promotor_regions(
    fasta_record: &fasta::Record,
    records: &[gtf::GTFRecord],
    transcripts: Vec<&gtf::GTFRecord>,
) -> Vec<PromotorRegion> {
    let ids: HashSet<_> =
        HashSet::from_iter(transcripts.iter().map(|&r| &r.attributes.transcript_id));

    let start_codons: Vec<_> = records
        .iter()
        .filter(|&r| r.feature_type == gtf::FeatureType::StartCodon)
        .filter(|&r| ids.contains(&r.attributes.transcript_id))
        .collect();

    let mut promotors: Vec<PromotorRegion> = Vec::with_capacity(start_codons.len());

    // ATG is at 1000, CAT at 98
    // in fasta: ATG at p.start+1000, CAT at p.end+98
    for start_codon in &start_codons {
        if start_codon.strand == gtf::Strand::Plus {
            let start = start_codon.start - 1001;
            let end = start_codon.start + 100;

            promotors.push(PromotorRegion {
                sequence: fasta_record.seq()[start..end].to_vec(),
                location: start,
                strand: gtf::Strand::Plus,
            });
        } else {
            let start = start_codon.end + 1000;
            let end = start_codon.end - 101;

            promotors.push(PromotorRegion {
                sequence: fasta_record.seq()[end..start].to_vec(),
                location: end,
                strand: gtf::Strand::Minus,
            });
        }
    }

    assert!(
        promotors.iter().map(|p| p.sequence.len()).max()
            == promotors.iter().map(|p| p.sequence.len()).min()
    );

    promotors
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
