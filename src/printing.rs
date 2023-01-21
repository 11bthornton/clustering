use std::collections::HashMap;

use itertools::Itertools;

use crate::helper_types::clustering_types::*;

pub fn print_umi_to_cdr(umi_to_cdr: &UMItoCDRClustering) -> HashMap<usize, usize> {

    let mut err_sizes_map = HashMap::new();

    let mut sorted = umi_to_cdr
        .iter()
        .collect_vec();

    sorted.sort_by(|one, two| {
        two.1
            .values()
            .max()
            .unwrap()
            .cmp(one.1.values().max().unwrap())
    });

    sorted.iter().for_each(|(key, values)| {

        let mut as_vec: Vec<_> = values.iter().collect_vec();
        as_vec.sort_by(|one, two| two.1.cmp(&one.1));

         // Will definitely have at least one...
        let (reference, ref_count) = as_vec[0];
        let copies = values.values().sum::<usize>();

        let reference_protein = reference.to_protein();

        println!("---- Clusters with identical barcode of {key} ({copies}) {{\n");
        println!("\tMajority   {reference}\t({ref_count} copies)\t{reference_protein}");

        for (presumed_error, count) in &as_vec[1..] {

            *err_sizes_map
                .entry(**count)
                .or_default() += 1;

            print!("\t           ");

            let _ = reference.pretty_align(presumed_error);
            let err_protein = presumed_error.to_protein();

            print!("\t({count} copies)\t");

            let _ = reference_protein.pretty_align(&err_protein);

            println!();
        }

        println!("\n}}\n");
    });

    err_sizes_map
}

pub fn print_cdr_to_umi(cdr_to_umi: &CDRtoUMIClustering, umi_to_cdr_opt: Option<Box<UMItoCDRClustering>>) -> Vec<f32> {
    // let mut err_sizes_map = HashMap::new();
    let mut hist_data = Vec::new();

    let mut sorted = cdr_to_umi
        .iter()
        .collect_vec();

    sorted.sort_by(|one, two| {
        two.1
            .values()
            .max()
            .unwrap()
            .cmp(one.1.values().max().unwrap())
    });

    println!("Total Count, Top Cluster Percentage, Second Cluster Percentage");

    sorted.iter().for_each(|(key, values)| {

        let mut as_vec: Vec<_> = values.iter().collect_vec();
        as_vec.sort_by(|one, two| two.1.cmp(&one.1));

         // Will definitely have at least one...
        let (reference, ref_count) = as_vec[0];
        let copies = values.values().sum::<usize>();

        let percentage = (*ref_count as f32 / copies as f32) * 100f32;

        // if percentage > 80 as f32 || (percentage < 10 as f32) {
        //     return;
        // };

        // println!("---- Clusters with identical cdr of {key} ({copies}) {{\n");
        // println!("\tMajority   {reference}   \t({ref_count} copies)\t{}%", percentage);

        // let sum_of_observations = values.values().sum::<usize>() as f32;
        // if copies > 20 {
        //     let percentage = (*ref_count as f32 / sum_of_observations) * 100f32;
            
        //     if as_vec.len() > 1 {
        //         let percentage2 = (*as_vec[1].1 as f32 / sum_of_observations) * 100f32;
        //         println!("{sum_of_observations}, {percentage}, {percentage2}");
        //     } else {
        //         println!("{sum_of_observations}, {percentage}, 0");
        //     }

            
        // }

        for (presumed_error, count) in &as_vec[1..] {

            // *err_sizes_map
            //     .entry(**count)
            //     .or_default() += 1;

            print!("\t           ");

            let _ = reference.normal_align(presumed_error);

            let percentage = (**count as f32 / values.values().sum::<usize>() as f32) * 100f32;
        

            print!("   \t({count} copies)\t{}%", percentage);

            // Print some extra information regarding the associated UMIs with this sequence.
            if let Some(umi_to_cdr) = &umi_to_cdr_opt {
                
                // -1 as know this umi is seen with the majority sequence.
                let other_unique_seq_count = &umi_to_cdr[*presumed_error].len() - 1;
                let count_seen_with = &umi_to_cdr[*presumed_error].values().sum::<usize>() - **count;

                let majority_count = &umi_to_cdr[*presumed_error]
                    .values()
                    .sorted()
                    .rev().nth(0).unwrap();

                let majority_sequence = *&umi_to_cdr[*presumed_error]
                    .iter()
                    .sorted_by(|one, two| two.1.cmp(&one.1))
                    .nth(0).unwrap().0;
                
                if other_unique_seq_count > 0 {
                    print!("\t(Seen with {other_unique_seq_count} other unique sequences, total of {count_seen_with} times)   \t");

                    if **count >= **majority_count {

                        print!("This UMI is most commonly associated with above CDR.");

                    } else {
                        print!(
                            "Max sequence cluster associated with this UMI is has size {} (/{count_seen_with})\t",
                            majority_count, 
                        );

                        key.pretty_align(majority_sequence);
                    }
                } else {
                    print!("\t(Only seen with this sequence)\t");
                }
            }

            println!();
        }

        println!("\n}}\n");
    });

    hist_data

}