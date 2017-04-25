# smmart-g2p

A prototype of the Genotype to phenotype user interface exists [here](https://g2p-ohsu.ddns.net/).   

![image](https://cloud.githubusercontent.com/assets/47808/24819664/66139700-1b9a-11e7-85ee-d49841486f0a.png)

## What is it?  Why use it?

* For researchers, who need to investigate genotype phenotype associations, smmart-g2p is a search tool that aggregates evidence from several knowledge bases unlike ad-hoc searches, the product allows the researcher to focus on the evidence, not on the search. [more](docs/smmart.pdf)

* Quickly determine the diseases, drugs and outcomes based on evidence from trusted sources. Find relevant articles and (soon) drug response data.

* Inform GA4GH G2P discussions


##  Where does the data come from?

Now:

* Jackson Lab [Clinical Knowledge Base](https://www.jax.org/clinical-genomics/clinical-offerings/ckb)
* Washington University [CIViC](https://civic.genome.wustl.edu/#/home)
* oncokb [Precision Oncology Knowledge Base](http://oncokb.org/#/)
* Cancer Genome Interpreter Cancer [bioMarkers database](https://www.cancergenomeinterpreter.org/biomarkers)
* GA4GH [reference server](https://github.com/ga4gh/ga4gh-server)

In  progress:

* [Molecular match](https://www.molecularmatch.com/technology.html#api-documentation)
* [BMEG](http://bmeg.compbio.ohsu.edu/)


## How to use it?

JUST GOOGLE IT:

* Use the search box like a google search. To search your data, enter your search criteria in the Query bar and press Enter or click Search to submit the request. For a full explanation of the search capabilities see [here](https://lucene.apache.org/core/2_9_4/queryparsersyntax.html)

* The charts and list are all tied to the search. Click to constrain your results


## Why are there a limited number of genes?

* We've constrained the data sets to a single use case from SMMART


## How do I import new data into it?

```
$ cd harvester
$ python harvester.py  -h
usage: harvester.py [-h] [--elastic_search ELASTIC_SEARCH]
                    [--elastic_index ELASTIC_INDEX] [--delete_index]
                    [--delete_source]
                    [--harvesters HARVESTERS [HARVESTERS ...]]

optional arguments:
  -h, --help            show this help message and exit
  --elastic_search ELASTIC_SEARCH, -es ELASTIC_SEARCH
                        elastic search endpoint
  --elastic_index ELASTIC_INDEX, -i ELASTIC_INDEX
                        elastic search index
  --delete_index, -d    delete elastic search index
  --delete_source, -ds  delete all content for source before harvest
  --harvesters HARVESTERS [HARVESTERS ...]
                        harvest from these sources. default: ['cgi_biomarkers', 'jax', 'civic', 'oncokb', 'g2p']
```

## How do I write a new harvester?
A `harvester` is a python module that implements this [duck typing](https://en.wikipedia.org/wiki/Duck_typing) interface.

```
#!/usr/bin/python


def harvest(genes):
    """ given a list of genes, yield an evidence item """
    # for gene in genes:
    #   gene_data = your_implementation_goes_here
    #      yield gene_data
    pass


def convert(gene_data):
    """ given a gene_data in it's original form, produce a feature_association """
    # gene: a string gene name
    # feature: a dict representing a ga4gh feature https://github.com/ga4gh/ga4gh-schemas/blob/master/src/main/proto/ga4gh/sequence_annotations.proto#L30
    # association: a dict representing a ga4gh g2p association https://github.com/ga4gh/ga4gh-schemas/blob/master/src/main/proto/ga4gh/genotype_phenotype.proto#L124
    #
    # feature_association = {'gene': gene ,
    #                        'feature': feature,
    #                        'association': association,
    #                        'source': 'my_source',
    #                        'my_source': {... original data from source ... }
    # yield feature_association
    pass


def harvest_and_convert(genes):
    """ get data from your source, convert it to ga4gh and return via yield """
    for gene_data in harvest(genes):
        for feature_association in convert(gene_data):
            yield feature_association

```

## How do I test it?

```
$ cd harvester
$ pytest -s -v
======================================================================================================================================================= test session starts ========================================================================================================================================================
platform darwin -- Python 2.7.13, pytest-3.0.7, py-1.4.33, pluggy-0.4.0 -- /usr/local/opt/python/bin/python2.7
cachedir: ../../.cache
rootdir: /Users/walsbr, inifile:
collected 13 items

tests/integration/test_elastic_silo.py::test_args PASSED
tests/integration/test_elastic_silo.py::test_init PASSED
tests/integration/test_elastic_silo.py::test_save PASSED
tests/integration/test_elastic_silo.py::test_delete_all PASSED
tests/integration/test_elastic_silo.py::test_delete_source PASSED
tests/integration/test_kafka_silo.py::test_populate_args PASSED
tests/integration/test_kafka_silo.py::test_init PASSED
tests/integration/test_kafka_silo.py::test_save PASSED
tests/integration/test_pb_deserialize.py::test_civic_pb PASSED
tests/integration/test_pb_deserialize.py::test_jax_pb PASSED
tests/integration/test_pb_deserialize.py::test_oncokb_pb PASSED
tests/integration/test_pb_deserialize.py::test_molecular_match_pb PASSED
tests/integration/test_pb_deserialize.py::test_cgi_pb PASSED
```

## What else do I need to know?

* See the README.md in harvester/tests/integration to see how harvested evidence is mapped to [protocol buffer messages](https://github.com/ohsu-comp-bio/bioschemas/blob/master/bioschemas/snapshot/proto/ohsu/g2p.proto).

## OK, I get it. But what about .... ?

### `NEXT STEPS`

* Work with users, gather feedback
* Load alternative data sources [literome, ensemble]
* Load smmart drugs [Olaparib, Folfox, Pembrolizumab, …]
* Integrate with bmeg (machine learning evidence)
* Improve data normalization
  * Variant naming (HGVS)
  * Ontologies (diseases, drugs, variants)
* Add GA4GH::G2P api  (or successor)
* Harden prototype:
  * python notebook
  * web app (deprecate kibana UI)


## Setup

* Create a .env file

```
ELASTIC_PORT=9200
KIBANA_PORT=5601
```

* update `services/nginx/default`  and `docker-compose.yml` to your certificate paths

```
# docker-compose.yml

      - "./compbio-tls:/compbio-tls"

# services/nginx/default

  ssl_certificate                 /compbio-tls/compbio_ohsu_edu_cert.cer;
  ssl_certificate_key             /compbio-tls/wild.compbio.ohsu.edu.key;
```

* create  `services/nginx/.htpasswd` 

Set to your userid:passwd.  See [here for an example](https://www.digitalocean.com/community/tutorials/how-to-set-up-basic-http-authentication-with-nginx-on-ubuntu-14-04#step-2-—-setting-up-http-basic-authentication-credentials) 

* load data
```
$ util/elastic-setup.sh
$ cd harvester; python harvester.py
```

* setup kibana

![image](https://cloud.githubusercontent.com/assets/47808/25396770/ace3daf6-299a-11e7-9300-07c6c885def6.png)



