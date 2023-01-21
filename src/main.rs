#![allow(incomplete_features)]
#![feature(array_chunks, generic_const_exprs)]

// use crate::affinity_csv::{find_set_difference, load_aff_data, threshold_experiment};

use crate::helper_types::{
    clustering_types::{CDRtoUMIClustering, UMItoCDRClustering},
    sequence::*,
};

use crate::printing::{print_cdr_to_umi, print_umi_to_cdr};

use std::{collections::HashMap, time::Instant};

use clap::Parser;
use seq_io::fastq::{Reader, Record};

use clargs::Config;

// mod affinity_csv;
mod clargs;
mod helper_types;
mod printing;

fn main() -> Result<(), ()> {

    let Config {
        print_type,
        take,
        qual
    } = Config::parse();

    println!("{:?}", Config::parse());

    // Collect arguments - num of pairs to process.

    // let affinity_data = load_aff_data("./data.csv")
    //     .expect("Problem with the csvs..");
    
    eprintln!("Extracted the bar codes from the csv.");

    let start = Instant::now();

    // Read from the files.
    let umi_reader = Reader::from_path("forward.fastq")
        .expect("Error with the UMIs.");

    let cdr_reader = Reader::from_path("reverse.fastq")
        .expect("Error with the CDRs.");

    // Zip the two iterators, coerce them into this library's expected types,
    // then accumulate them into two hashmaps.
    //
    // This will create two hashmaps in one:
    // - The first being a clustering of unique UMI -> CDR.
    // - The second being a clustering of unique CDR -> UMI.
    let mut count = 0;
    let (mut umi_to_cdr, cdr_to_umi): (UMItoCDRClustering, CDRtoUMIClustering) = umi_reader
        .into_records()
        .take(take as usize)
        .zip(cdr_reader.into_records().take(take as usize))
        .map(|(umi, cdr)| (UMI::try_from((umi.as_ref().unwrap().seq(), umi.as_ref().unwrap().qual())).unwrap(), CDR::try_from((cdr.as_ref().unwrap().seq(), cdr.as_ref().unwrap().qual())).unwrap()))
        .fold(
            (HashMap::new().into(), HashMap::new().into()),
            |(mut umicdrmap, mut cdrumimap), (umi, cdr)| {
                count += 1;

                if count % 1000000 == 0 {
                    eprintln!(
                        "{:03} million reads processed in {:?}",
                        count / 1000000,
                        start.elapsed()
                    );
                }

                *umicdrmap
                    .entry(umi.clone())
                    .or_default()
                    .entry(cdr.clone())
                    .or_insert(0) += if cdr.ratio_ns() < 0.55 { 1 } else { 0 };

                *cdrumimap
                    .entry(cdr)
                    .or_default()
                    .entry(umi)
                    .or_insert(0) += 1;

                (umicdrmap, cdrumimap)
            },
        );

    eprintln!(
        "Completed initial clustering of {count} total pairs (both ways) on identical {print_type} in {:?}. {} total umi clusters.",
        start.elapsed(),
        umi_to_cdr.len()
    );

    // let exp_results = threshold_experiment(&affinity_data, &mut umi_to_cdr, 2500);

    // println!("[");
    // exp_results.iter()
    //     .for_each(|(threshold, (count, missing))| {
    //         // println!("When keep-threshold is > {threshold}, {count} umis from the affinity data are not found in the index {:?}", missing);
            
    //         println!("({},{}),", threshold, count)
    //   });
    // println!("]");

    // eprintln!("Before filtering {difference_before} were not found, but after {difference_after} were not found (out of {})", affinity_data.len());

    if print_type.as_str() == "cdr" {
        print_umi_to_cdr(&umi_to_cdr);
    } else {

        let res = print_cdr_to_umi(&cdr_to_umi, Some(Box::new(umi_to_cdr)));
        // let res = print_cdr_to_umi(&cdr_to_umi, None);

        // println!("{:?}", res);

    }

    Ok(())
}