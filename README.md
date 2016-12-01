# Installation

```{r}
devtools::install_github('seandavi/VCFWrenchR')
```

# Very basic usage

Read in a VCF file using VariantAnnotation.

```{r}
library(VariantAnnotation)
vcffile = system.file(package='VariantAnnotation',path='extdata/chr22.vcf.gz')
v = readVcf(vcffile,'hg19')
head(v)
````

## Convert to data.frame

```{r}
head(as.data.frame(v))
```

## Convert to JSON

```{r}
as.json(head(v,2),pretty=TRUE)
```

## Write line-delimited json

```{r}
con = tempfile()
jsonlite::stream_out(as.data.frame(v),con)
close(con)
```

# Working with `jq` and ElasticSearch

## Load line-delimited json variants into ElasticSearch

TODO: Do this in R....

```{sh}
# assumes jq installed
cat abc.json | \
    jq -c '. | {"index": {"_index": "variants", "_type": "variant"}}, .' | \
    curl -XPOST localhost:9200/_bulk --data-binary @-
```
