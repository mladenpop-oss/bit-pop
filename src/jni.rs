use jni::objects::{JClass, JString};
use jni::sys::{jint, jstring};
use jni::JNIEnv;

use crate::{BitPop, fastq};

#[no_mangle]
pub extern "system" fn Java_com_bitpop_MainActivity_mapReads(
    mut env: JNIEnv,
    _: JClass,
    index_path: JString,
    reads_path: JString,
    output_path: JString,
) -> jstring {
    let result = (|| -> Result<String, Box<dyn std::error::Error>> {
        let idx: String = env.get_string(&index_path)?.into();
        let rds: String = env.get_string(&reads_path)?.into();
        let out: String = env.get_string(&output_path)?.into();

        let bp = BitPop::deserialize_from_file(&idx)?;
        let reads = fastq::parse_fastq(&rds)?;

        let mut results = Vec::new();
        for (name, seq, _) in reads {
            let mapped = bp.map_read(&seq, 10);
            if let Some(best) = mapped.first() {
                let genome_name = bp.genome_name(best.genome_id).unwrap_or("unknown");
                results.push(format!(
                    "{}\t{}\t{:.2}",
                    name, genome_name, best.score
                ));
            }
        }

        std::fs::write(&out, results.join("\n") + "\n")?;
        Ok(format!("Mapped {} reads to {}", results.len(), out))
    })();

    let msg = match result {
        Ok(m) => m,
        Err(e) => format!("Error: {}", e),
    };

    env.new_string(&msg)
        .unwrap_or_else(|_| env.new_string("Unknown error").unwrap())
        .into_raw()
}

#[no_mangle]
pub extern "system" fn Java_com_bitpop_MainActivity_buildIndex(
    mut env: JNIEnv,
    _: JClass,
    fasta_path: JString,
    output_path: JString,
    k_mer_size: jint,
) -> jstring {
    let result = (|| -> Result<String, Box<dyn std::error::Error>> {
        let fasta: String = env.get_string(&fasta_path)?.into();
        let output: String = env.get_string(&output_path)?.into();

        let mut bp = BitPop::new(k_mer_size as usize);
        let genomes = fastq::parse_fasta(&fasta)?;

        for (name, seq) in genomes {
            bp.add_genome(&name, &seq);
        }

        bp.build();
        bp.serialize_to_file(&output)?;

        Ok(format!("Index saved to {}", output))
    })();

    let msg = match result {
        Ok(m) => m,
        Err(e) => format!("Error: {}", e),
    };

    env.new_string(&msg)
        .unwrap_or_else(|_| env.new_string("Unknown error").unwrap())
        .into_raw()
}

#[no_mangle]
pub extern "system" fn Java_com_bitpop_MainActivity_getGenomeNames(
    mut env: JNIEnv,
    _: JClass,
    index_path: JString,
) -> jstring {
    let result = (|| -> Result<String, Box<dyn std::error::Error>> {
        let idx: String = env.get_string(&index_path)?.into();
        let bp = BitPop::deserialize_from_file(&idx)?;
        let names: Vec<String> = bp.genome_names_ordered();
        Ok(names.join(","))
    })();

    let msg = match result {
        Ok(m) => m,
        Err(e) => format!("Error: {}", e),
    };

    env.new_string(&msg)
        .unwrap_or_else(|_| env.new_string("Unknown error").unwrap())
        .into_raw()
}
