import vcf

def filter_vcf(input_vcf, output_vcf):
    """Filter VCF based on certain criteria (e.g., quality)"""
    vcf_reader = vcf.Reader(open(input_vcf, 'r'))
    vcf_writer = vcf.Writer(open(output_vcf, 'w'), vcf_reader)

    for record in vcf_reader:
        if record.QUAL >= 30.0:  # Example: Filter variants with quality >= 30
            vcf_writer.write_record(record)
    print(f"Filtered VCF written to {output_vcf}")

# Example usage:
if __name__ == "__main__":
    input_vcf = "/path/to/input.vcf"
    output_vcf = "/path/to/output_filtered.vcf"
    filter_vcf(input_vcf, output_vcf)
