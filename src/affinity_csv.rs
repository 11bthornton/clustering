use std::collections::{HashSet, HashMap};

use csv;
use itertools::Itertools;

use crate::helper_types::{sequence::UMI, clustering_types::UMItoCDRClustering};



pub fn load_aff_data(path: &str) -> Result<HashSet<UMI>, csv::Error>  {
    let mut data = csv::Reader::from_path(path)?;

    Ok(
        data
        .byte_records()
        .map(|record| UMI::try_from(&record.unwrap()[1]).expect("Barcode not of correct type"))
        .fold(HashSet::new(), |mut set, record| {

            set.insert(record);

            set
        })
    )
}

pub fn find_set_difference(affinity_umis: &HashSet<UMI>, umi_index: &UMItoCDRClustering) -> (usize, HashSet<UMI>) {

    let mut count_not_found = 0;
    let mut not_found = HashSet::new();

    affinity_umis
        .iter()
        .for_each(|umi| {
            if !umi_index.contains_key(umi) {
                // println!("{umi} not found in index");
                count_not_found += 1;

                not_found.insert(umi.clone());
            }
        });
    
    (count_not_found, not_found)
}

pub fn threshold_experiment(affinity_umis: &HashSet<UMI>, umi_index: &mut UMItoCDRClustering, threshold: usize) -> Vec<(usize, (usize, HashSet<UMI>))> {

    let mut hashmap = HashMap::new();

    for threshold in 0..(threshold + 1) {

        eprintln!("doing {}", threshold);

        umi_index.filter_threshold(threshold);

        hashmap.insert(
            threshold, 
            find_set_difference(affinity_umis, umi_index)
        );
    }

    let mut vec = hashmap
        .into_iter()
        .collect_vec();

    vec.sort_by(|(ind1, _), (ind2, _)| ind1.cmp(ind2) );


    vec
}