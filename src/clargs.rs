use clap::Parser;

#[derive(Debug, Parser)]
pub struct Config {

    #[arg(short, long, default_value_t = String::from("umi"))]
    pub print_type: String,

    #[arg(short, long, default_value_t = 180000000)]
    pub take: u64,

    #[arg(short, long, default_value_t = false)]
    pub qual: bool

}