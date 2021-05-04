#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    "nf-core/gemmaker": ["v_pipeline.txt", r"(\S+)"],
    "Nextflow": ["v_nextflow.txt", r"(\S+)"],
    "FastQC": ["v_fastqc.txt", r"FastQC v(\S+)"],
    "MultiQC": ["v_multiqc.txt", r"multiqc, version (\S+)"],
    "Kallisto": ["v_kallisto.txt", r"kallisto, version (\S+)"],
    "Hisast2": ["v_hisat2.txt", r".+ version (\S+)"],
    "Salmon": ["v_salmon.txt", r"salmon (\S+)"],
    "Python": ["v_python.txt", r"Python (\S+)"],
    "SAMTools": ["v_samtools.txt", r"samtools (\S+)"],
    "fastq-dump": ["v_fastq_dump.txt", r"\"fastq-dump\" version (\S+)"],
    "Stringtie": ["v_stringtie.txt", r"(\S+)"],
    "Trimmomatic": ["v_trimmomatic.txt", r"(\S+)"],
    "Pandas": ["v_pandasc.txt", r"pandas==(\S+)"]
}
results = OrderedDict()
results["nf-core/gemmaker"] = '<span style="color:#999999;">N/A</span>'
results["Nextflow"] = '<span style="color:#999999;">N/A</span>'
results["FastQC"] = '<span style="color:#999999;">N/A</span>'
results["MultiQC"] = '<span style="color:#999999;">N/A</span>'
results["Kallisto"] = '<span style="color:#999999;">N/A</span>'
results["Hisast2"] = '<span style="color:#999999;">N/A</span>'
results["Salmon"] = '<span style="color:#999999;">N/A</span>'
results["Python"] = '<span style="color:#999999;">N/A</span>'
results["SAMTools"] = '<span style="color:#999999;">N/A</span>'
results["fastq-dump"] = '<span style="color:#999999;">N/A</span>'
results["Stringtie"] = '<span style="color:#999999;">N/A</span>'
results["Trimmomatic"] = '<span style="color:#999999;">N/A</span>'
results["Pandas"] = '<span style="color:#999999;">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del results[k]

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'nf-core/gemmaker Software Versions'
section_href: 'https://github.com/nf-core/gemmaker'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
for k, v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.csv", "w") as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))
