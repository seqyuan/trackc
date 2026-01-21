class BedBuilder:
    __slots__ = (
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "thick_start",
        "thick_end",
        "rgb",
        "block_count",
        "exons",
        "biotype",
    )

    def __init__(self, chrom, start, end, name, strand, biotype):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.score = 0
        self.strand = strand
        self.thick_start = "."
        self.thick_end = "."
        self.rgb = 0
        self.block_count = 0
        self.exons = []
        self.biotype = biotype

    def add_exon(self, start, end):
        self.exons.append((int(start), int(end)))

    def to_string(self, include_biotype):
        if not self.exons:
            return ""

        block_count = len(self.exons)
        self.exons.sort()
        block_sizes = [str(e[1] - e[0]) for e in self.exons]
        block_starts = [str(e[0]) for e in self.exons]

        cols = [
            self.chrom,
            str(self.start),
            str(self.end),
            self.name,
            str(self.score),
            self.strand,
            str(self.thick_start),
            str(self.thick_end),
            str(self.rgb),
            str(block_count),
            ",".join(block_sizes),
            ",".join(block_starts),
        ]
        if include_biotype:
            cols.append(self.biotype)

        return "\t".join(cols)


def parse_attributes(attr_str, sep=";"):
    attrs = {}
    for item in attr_str.split(sep):
        if not item:
            continue
        parts = item.strip().split(" ", 1)
        if len(parts) == 2:
            attrs[parts[0]] = parts[1].strip('"')
    return attrs


def gtf2bed(
    gtf_file,
    outfile,
    sep=";",
    gene_name_tag="gene_name",
    gene_id_tag="gene_id",
    gene_biotype_tag="gene_biotype",
    biotype2bed13=False,
):
    current_gene = None
    with open(gtf_file, "r") as f_in, open(outfile, "w") as f_out:
        for line in f_in:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            feature_type = parts[2]

            if feature_type == "gene":
                if current_gene:
                    f_out.write(current_gene.to_string(biotype2bed13) + "\n")

                attrs = parse_attributes(parts[8], sep)

                name = attrs.get(gene_name_tag) or attrs.get(gene_id_tag) or "."
                if name == ".":
                    print(f"Warning: No name found for line: {parts[8]}")

                biotype = attrs.get(gene_biotype_tag) or attrs.get("gene_type") or "."

                current_gene = BedBuilder(
                    chrom=parts[0],
                    start=parts[3],
                    end=parts[4],
                    name=name,
                    strand=parts[6],
                    biotype=biotype,
                )

            elif feature_type == "exon":
                if current_gene:
                    current_gene.add_exon(parts[3], parts[4])

        if current_gene:
            f_out.write(current_gene.to_string(biotype2bed13) + "\n")

    print("gtf2bed finished")
