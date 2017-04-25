__author__ = 'andra'

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../ProteinBoxBot_Core")
from pprint import pprint
import requests
from SPARQLWrapper import SPARQLWrapper, JSON
import ProteinBoxBot_Core.PBB_Core as PBB_Core

import ProteinBoxBot_Core.PBB_login as PBB_login
from time import gmtime, strftime
import copy
import traceback
import time

logincreds = PBB_login.WDLogin("ProteinBoxBot", os.environ['wikidataApi'])

# chromosomes
#Chromosomes dict
chromosomes = dict()
chromosomes['1'] = "Q430258"
chromosomes['2'] = "Q638893"
chromosomes['3'] = "Q668633"
chromosomes['4'] = "Q836605"
chromosomes['5'] = "Q840741"
chromosomes['6'] = "Q540857"
chromosomes['7'] = "Q657319"
chromosomes['8'] = "Q572848"
chromosomes['9'] = "Q840604"
chromosomes['10'] = "Q840737"
chromosomes['11'] = "Q847096"
chromosomes['12'] = "Q847102"
chromosomes['13'] = "Q840734"
chromosomes['14'] = "Q138955"
chromosomes['15'] = "Q765245"
chromosomes['16'] = "Q742870"
chromosomes['17'] = "Q220677"
chromosomes['18'] = "Q780468"
chromosomes['19'] = "Q510786"
chromosomes['20'] = "Q666752"
chromosomes['21'] = "Q753218"
chromosomes['22'] = "Q753805"
chromosomes['22'] = "Q753805"
chromosomes['X'] = "Q61333"
chromosomes['Y'] = "Q202771"
chromosomes["MT"] = "Q27075"

r = requests.get('https://civic.genome.wustl.edu/api/variants?count=10000')
variant_data = r.json()


ignore_synonym_list = [
    "AMPLIFICATION",
    "EXPRESSION",
    "DELETION",
    "LOSS",
    "LOSS-OF-FUNCTION",
    "MUTATION",
    "NUCLEAR EXPRESSION",
    "OVEREXPRESSION",
    "UNDEREXPRESSION",
    "3\' UTR MUTATION",
    "BIALLELIC INACTIVATION",
    "EWSR1-FLI1",
    "EXON 12 MUTATION",
    "EXON 9 MUTATION",
    "FRAMESHIFT TRUNCATION",
    "G12",
    "G12/G13",
    "G13D",
    "METHYLATION",
    "PHOSPHORYLATION",
    "PROMOTER HYPERMETHYLATION",
    "PROMOTER METHYLATION",
    "SERUM LEVELS",
    "TMPRSS2-ERG",
    "TRUNCATING MUTATION",
]

for record in variant_data['records']:
 # try:
    fast_run_base_filter = {'P3329': ''}
    fast_run = True
    print(record['id'])
    # Reference section
    # Prepare references
    # variant_id = "499"
    # variant_id = "183"
    variant_id = str(record['id'])
    # variant_id = "12"
    refStatedIn = PBB_Core.WDItemID(value="Q27612411", prop_nr='P248', is_reference=True)
    timeStringNow = strftime("+%Y-%m-%dT00:00:00Z", gmtime())
    refRetrieved = PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)
    refReferenceURL = PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/variants/"+variant_id, prop_nr="P854", is_reference=True)
    variant_reference = [refStatedIn, refRetrieved, refReferenceURL]

    genomeBuildQualifier = PBB_Core.WDItemID(value="Q21067546", prop_nr='P659', is_qualifier=True)


    r = requests.get('https://civic.genome.wustl.edu/api/variants/'+variant_id)
    variant_data = r.json()
    pprint(variant_data)

    prep = dict()

    sparql = SPARQLWrapper("https://query.wikidata.org/bigdata/namespace/wdq/sparql")

    ncbi_geneQuery = """SELECT ?item ?itemLabel
        WHERE
        {
            ?item wdt:P351 """
    ncbi_geneQuery += "\""+str(variant_data["entrez_id"])
    ncbi_geneQuery += """\" .
            SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
        }
        """
    print(ncbi_geneQuery)
    sparql.setQuery(ncbi_geneQuery)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    for result in results["results"]["bindings"]:
       prep['P3433'] = [PBB_Core.WDItemID(value=result["item"]["value"].replace("http://www.wikidata.org/entity/", ""), prop_nr='P3433', references=[copy.deepcopy(variant_reference)])]
       print("wd entrez:" + result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))


    # part of
    partof = variant_data["entrez_id"]

    # variant_id
    prep['P3329'] = [PBB_Core.WDString(value=variant_id, prop_nr='P3329', references=[copy.deepcopy(variant_reference)])]


    #coordinates
    coordinates = variant_data["coordinates"]
    if coordinates["chromosome"] != None:
        prep['P1057'] = [PBB_Core.WDItemID(value=chromosomes[coordinates["chromosome"]], prop_nr='P1057', references=[copy.deepcopy(variant_reference)], qualifiers=[copy.deepcopy(genomeBuildQualifier)])]
        if coordinates["chromosome2"] != None:
            prep['P1057'].append(PBB_Core.WDItemID(value=chromosomes[coordinates["chromosome2"]], prop_nr='P1057', references=[copy.deepcopy(variant_reference)], qualifiers=[copy.deepcopy(genomeBuildQualifier)]))

        # genomic start
        prep['P644'] = [PBB_Core.WDString(value=str(coordinates["start"]), prop_nr='P644', references=[copy.deepcopy(variant_reference)], qualifiers=[copy.deepcopy(genomeBuildQualifier)])]
        prep['P645'] = [PBB_Core.WDString(value=str(coordinates["stop"]), prop_nr='P645', references=[copy.deepcopy(variant_reference)], qualifiers=[copy.deepcopy(genomeBuildQualifier)])]

        if coordinates["start2"] != None:
            prep['P644'].append(PBB_Core.WDString(value=str(coordinates["start2"]), prop_nr='P644', references=[copy.deepcopy(variant_reference)], qualifiers=[copy.deepcopy(genomeBuildQualifier)]))
            prep['P645'].append(PBB_Core.WDString(value=str(coordinates["stop2"]), prop_nr='P645', references=[copy.deepcopy(variant_reference)], qualifiers=[copy.deepcopy(genomeBuildQualifier)]))

    query = """
            SELECT DISTINCT ?item  ?itemLabel ?alias
            WHERE
            {
                ?item p:P528 ?o .
                ?o pq:P972	wd:Q7452458
                OPTIONAL {?item skos:altLabel ?alias .
                		  FILTER (lang(?alias)="en")}

                SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
              }
            """
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    seqO = dict()
    for result in results["results"]["bindings"]:
        seqO[result["itemLabel"]["value"]] = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
        if "alias" in result.keys():
          if result["alias"]["value"] != "":
            seqO[result["alias"]["value"]] = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
        seqO[result["itemLabel"]["value"]] = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
        print("wd disease:" + result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))

    prep["P31"] = []
    for variant_type in variant_data["variant_types"]:
        if variant_type["name"] == "N/A":
            continue
        prep['P31'].append(PBB_Core.WDItemID(value=seqO[variant_type["display_name"]], prop_nr='P31', references=[copy.deepcopy(variant_reference)]))

        print(variant_type["display_name"])

    name = variant_data["variant_aliases"]
    hgvs_expressions = variant_data["hgvs_expressions"]

    print(name)

    evidence = dict()
    evidence["P3354"]=dict()
    evidence["P3355"]=dict()
    evidence["P3356"]=dict()
    evidence["P3357"]=dict()
    evidence["P3358"]=dict()
    evidence["P3359"]=dict()
    for evidence_item in variant_data["evidence_items"]:
        #if evidence_item["evidence_level"] == "A" or evidence_item["evidence_level"] == "B" or evidence_item["evidence_level"] == "C" :
            ## Disease
            Disease = None
            if evidence_item["disease"]["doid"] != None or evidence_item["disease"]["doid"] != "":
                print("DOID:"+evidence_item["disease"]["doid"])

                query = """SELECT ?item ?itemLabel
                    WHERE
                    {
                        ?item wdt:P699 """
                query += "\"DOID:"+evidence_item["disease"]["doid"]
                query += """\" .
                        SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
                    }
                    """
                sparql.setQuery(query)
                sparql.setReturnFormat(JSON)
                results = sparql.query().convert()
                for result in results["results"]["bindings"]:
                  print("wd disease:" + result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))
                  disease = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
            else:
                continue

            ## Drugs
            print(evidence_item["drugs"])
            wd_drugs = []
            for drug in evidence_item["drugs"]:
                query = """
                SELECT ?item ?itemLabel
                WHERE
                {
                    ?item wdt:P31 wd:Q12140 ;
                          rdfs:label ?label .
                    FILTER regex(?label,  \""""
                query += drug["name"]
                query += """\", \"i\")
                }
                """
                print(query)
                sparql.setQuery(query)
                sparql.setReturnFormat(JSON)
                results = sparql.query().convert()
                for result in results["results"]["bindings"]:
                    print("wd drug:" + result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))
                    wd_drugs.append(result["item"]["value"].replace("http://www.wikidata.org/entity/", ""))

            ## Pubmed
            print("PMID: "+evidence_item["source"]["pubmed_id"])
            pubmedQuery = """SELECT ?item ?itemLabel
                WHERE
                {
                    ?item wdt:P698 """
            pubmedQuery += "\""+evidence_item["source"]["pubmed_id"]
            pubmedQuery += """\" .
                    SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
                }
                """

            print(pubmedQuery)
            sparql.setQuery(pubmedQuery)
            sparql.setReturnFormat(JSON)
            results = sparql.query().convert()
            references = []

            for result in results["results"]["bindings"]:
               pubmed_entry = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
               refStatedIn = PBB_Core.WDItemID(value=pubmed_entry, prop_nr='P248', is_reference=True)
               references.append(refStatedIn)

            timeStringNow = strftime("+%Y-%m-%dT00:00:00Z", gmtime())
            refRetrieved = PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)
            #pubmed_references.append(refRetrieved)
            #refReferenceURL = PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/variants/"+variant_id, prop_nr="P854", is_reference=True)
            #pubmed_references.append(refReferenceURL)

            # Positive therapeutic predictor
            if evidence_item["evidence_type"] ==  "Predictive" and evidence_item["clinical_significance"] == "Sensitivity" and evidence_item["evidence_direction"] == "Supports":
                for drug in wd_drugs:
                    if drug not in evidence["P3354"].keys():
                        evidence["P3354"][drug] = dict()
                    if disease not in evidence["P3354"][drug].keys():
                        evidence["P3354"][drug][disease] = dict()
                    if "stated_in" not in evidence["P3354"][drug][disease].keys():
                        evidence["P3354"][drug][disease]["stated_in"] = dict()
                    if "references" not in evidence["P3354"][drug][disease]["stated_in"].keys():
                        evidence["P3354"][drug][disease]["stated_in"]["references"]=[]
                    if "id" not in evidence["P3354"][drug][disease]["stated_in"].keys():
                        evidence["P3354"][drug][disease]["stated_in"]["id"]=[]
                    if str(evidence_item["id"]) not in evidence["P3354"][drug][disease]["stated_in"]["id"]:
                        evidence["P3354"][drug][disease]["stated_in"]["id"].append(str(evidence_item["id"]))
                    if pubmed_entry not in evidence["P3354"][drug][disease]["stated_in"]["references"]:
                        evidence["P3354"][drug][disease]["stated_in"]["references"].append(pubmed_entry)

            # Positive therapeutic predictor disputed by
            if evidence_item["evidence_type"] ==  "Predictive" and evidence_item["clinical_significance"] == "Sensitivity" and evidence_item["evidence_direction"] == "Does Not Support":
                for drug in wd_drugs:
                    if drug not in evidence["P3354"].keys():
                        evidence["P3354"][drug] = dict()
                    if disease not in evidence["P3354"][drug].keys():
                        evidence["P3354"][drug][disease] = dict()
                    if "disputed_by" not in evidence["P3354"][drug][disease].keys():
                        evidence["P3354"][drug][disease]["disputed_by"] = dict()
                    if "references" not in evidence["P3354"][drug][disease]["disputed_by"].keys():
                        evidence["P3354"][drug][disease]["disputed_by"]["references"]=[]
                    if "id" not in evidence["P3354"][drug][disease]["disputed_by"].keys():
                        evidence["P3354"][drug][disease]["disputed_by"]["id"]=[]
                    if str(evidence_item["id"]) not in evidence["P3354"][drug][disease]["disputed_by"]["id"]:
                        evidence["P3354"][drug][disease]["disputed_by"]["id"].append(str(evidence_item["id"]))
                    if pubmed_entry not in evidence["P3354"][drug][disease]["disputed_by"]["references"]:
                        evidence["P3354"][drug][disease]["disputed_by"]["references"].append(pubmed_entry)

            # Negative therapeutic predictor
            if evidence_item["evidence_type"] ==  "Predictive" and evidence_item["clinical_significance"] == "Resistance or Non-Response" and evidence_item["evidence_direction"] == "Supports":
                for drug in wd_drugs:
                    if drug not in evidence["P3355"].keys():
                        evidence["P3355"][drug] = dict()
                    if disease not in evidence["P3355"][drug].keys():
                        evidence["P3355"][drug][disease] = dict()
                    if "stated_in" not in evidence["P3355"][drug][disease].keys():
                        evidence["P3355"][drug][disease]["stated_in"] = dict()
                    if "references" not in evidence["P3355"][drug][disease]["stated_in"].keys():
                        evidence["P3355"][drug][disease]["stated_in"]["references"]=[]
                    if "id" not in evidence["P3355"][drug][disease]["stated_in"].keys():
                        evidence["P3355"][drug][disease]["stated_in"]["id"]=[]
                    if str(evidence_item["id"]) not in evidence["P3355"][drug][disease]["stated_in"]["id"]:
                        evidence["P3355"][drug][disease]["stated_in"]["id"].append(str(evidence_item["id"]))
                    if pubmed_entry not in evidence["P3355"][drug][disease]["stated_in"]["references"]:
                        evidence["P3355"][drug][disease]["stated_in"]["references"].append(pubmed_entry)

            # Negative therapeutic predictor disputed by
            if evidence_item["evidence_type"] ==  "Predictive" and evidence_item["clinical_significance"] == "Resistance or Non-Response" and evidence_item["evidence_direction"] == "Does Not Support":
                for drug in wd_drugs:
                    if drug not in evidence["P3355"].keys():
                        evidence["P3355"][drug] = dict()
                    if disease not in evidence["P3355"][drug].keys():
                        evidence["P3355"][drug][disease] = dict()
                    if "disputed_by" not in evidence["P3355"][drug][disease].keys():
                        evidence["P3355"][drug][disease]["disputed_by"] = dict()
                    if "references" not in evidence["P3355"][drug][disease]["disputed_by"].keys():
                        evidence["P3355"][drug][disease]["disputed_by"]["references"]=[]
                    if "id" not in evidence["P3355"][drug][disease]["disputed_by"].keys():
                        evidence["P3355"][drug][disease]["disputed_by"]["id"]=[]
                    if str(evidence_item["id"]) not in  evidence["P3355"][drug][disease]["disputed_by"]["id"]:
                        evidence["P3355"][drug][disease]["disputed_by"]["id"].append(str(evidence_item["id"]))
                    if pubmed_entry not in evidence["P3355"][drug][disease]["disputed_by"]["references"]:
                        evidence["P3355"][drug][disease]["disputed_by"]["references"].append(pubmed_entry)

            # Positive diagnostic predictor (stated in)
            if evidence_item["evidence_type"] ==  "Diagnostic" and evidence_item["clinical_significance"] == "Positive" and evidence_item["evidence_direction"] == "Supports":
                if disease not in evidence["P3356"].keys():
                    evidence["P3356"][disease] = dict()
                if "references" not in evidence["P3356"][disease].keys():
                        evidence["P3356"][disease]["references"]=dict()
                if "stated_in" not in evidence["P3356"][disease]["references"].keys():
                        evidence["P3356"][disease]["references"]["stated_in"]=[]
                evidence["P3356"][disease]["id"] = str(evidence_item["id"])
                evidence["P3356"][disease]["references"]["stated_in"].append(pubmed_entry)

            # Positive diagnostic predictor (disputed by)
            if evidence_item["evidence_type"] ==  "Diagnostic" and evidence_item["clinical_significance"] == "Positive" and evidence_item["evidence_direction"] == "Does Not Support":
                if disease not in evidence["P3356"].keys():
                    evidence["P3356"][disease] = dict()
                if "references" not in evidence["P3356"][disease].keys():
                        evidence["P3356"][disease]["references"]=dict()
                if "disputed_by" not in evidence["P3356"][disease]["references"].keys():
                        evidence["P3356"][disease]["references"]["disputed_by"]=[]
                evidence["P3356"][disease]["id"] = str(evidence_item["id"])
                evidence["P3356"][disease]["references"]["disputed_by"].append(pubmed_entry)

            # Negative diagnostic predictor (stated in)
            if evidence_item["evidence_type"] ==  "Diagnostic" and evidence_item["clinical_significance"] == "Negative" and evidence_item["evidence_direction"] == "Supports":
                if disease not in evidence["P3357"].keys():
                    evidence["P3357"][disease] = dict()
                if "references" not in evidence["P3357"][disease].keys():
                        evidence["P3357"][disease]["references"]=dict()
                if "stated_in" not in evidence["P3357"][disease]["references"].keys():
                        evidence["P3357"][disease]["references"]["stated_in"]=[]
                evidence["P3357"][disease]["id"] = str(evidence_item["id"])
                evidence["P3357"][disease]["references"]["stated_in"].append(pubmed_entry)

            # Negative diagnostic predictor (disputed by)
            if evidence_item["evidence_type"] ==  "Diagnostic" and evidence_item["clinical_significance"] == "Negative" and evidence_item["evidence_direction"] == "Does Not Support":
                if disease not in evidence["P3357"].keys():
                    evidence["P3357"][disease] = dict()
                if "references" not in evidence["P3357"][disease].keys():
                        evidence["P3357"][disease]["references"]=dict()
                if "disputed_by" not in evidence["P3357"][disease]["references"].keys():
                        evidence["P3357"][disease]["references"]["disputed_by"]=[]
                evidence["P3357"][disease]["id"] = str(evidence_item["id"])
                evidence["P3357"][disease]["references"]["disputed_by"].append(pubmed_entry)

            # Positive prognostic predictor
            if evidence_item["evidence_type"] ==  "Prognostic" and evidence_item["clinical_significance"] == "Better Outcome" and evidence_item["evidence_direction"] == "Supports":
                if disease not in evidence["P3358"].keys():
                    evidence["P3358"][disease] = dict()
                if "references" not in evidence["P3358"][disease].keys():
                        evidence["P3358"][disease]["references"]=dict()
                if "stated_in" not in evidence["P3358"][disease]["references"].keys():
                        evidence["P3358"][disease]["references"]["stated_in"]=[]
                evidence["P3358"][disease]["id"] = str(evidence_item["id"])
                evidence["P3358"][disease]["references"]["stated_in"].append(pubmed_entry)


            # Negative prognostic predictor (disputed by)
            if evidence_item["evidence_type"] ==  "Prognostic" and evidence_item["clinical_significance"] == "Poor Outcome" and evidence_item["evidence_direction"] == "Does Not Support":
                if disease not in evidence["P3359"].keys():
                    evidence["P3359"][disease] = dict()
                    if "references" not in evidence["P3359"][disease].keys():
                        evidence["P3359"][disease]["references"]=dict()
                    if "disputed_by" not in evidence["P3359"][disease]["references"].keys():
                        evidence["P3359"][disease]["references"]["disputed_by"]=[]
                    evidence["P3359"][disease]["id"] = str(evidence_item["id"])
                    evidence["P3359"][disease]["references"]["disputed_by"].append(pubmed_entry)



    pprint(evidence)

    if len(evidence["P3354"]) > 0:
        prep["P3354"] = []
        for drug in evidence["P3354"].keys():
            for disease in evidence["P3354"][drug].keys():
                if disease != None:
                    references = [PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)]
                    disp_references = [PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)]
                    disease_qualifier = [PBB_Core.WDItemID(value=disease, prop_nr='P2175', is_qualifier=True)]
                    if "stated_in" in evidence["P3354"][drug][disease].keys():
                        for statedin in evidence["P3354"][drug][disease]["stated_in"]["references"]:
                            references.append(PBB_Core.WDItemID(value=statedin, prop_nr='P248', is_reference=True))
                        for evidence_id in evidence["P3354"][drug][disease]["stated_in"]["id"]:
                            references.append(PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/evidence/"+evidence_id, prop_nr="P854", is_reference=True))
                        prep["P3354"].append(PBB_Core.WDItemID(value=drug, prop_nr='P3354', references=[copy.deepcopy(references)] , qualifiers=copy.deepcopy(disease_qualifier)))

                    if "disputed_by" in evidence["P3354"][drug][disease].keys():
                        for disputedby in evidence["P3354"][drug][disease]["disputed_by"]["references"]:
                            disease_qualifier.append(PBB_Core.WDItemID(value=disputedby, prop_nr='P1310', is_qualifier=True))
                            disp_references.append(PBB_Core.WDItemID(value=disputedby, prop_nr='P248', is_reference=True))
                        for evidence_id in evidence["P3354"][drug][disease]["disputed_by"]["id"]:
                            disp_references.append(PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/evidence/"+evidence_id, prop_nr="P854", is_reference=True))
                        prep["P3354"].append(PBB_Core.WDItemID(value=drug, prop_nr='P3354', references=[copy.deepcopy(disp_references)] , qualifiers=copy.deepcopy(disease_qualifier)))

    if len(evidence["P3355"]) > 0:
        prep["P3355"] = []
        for drug in evidence["P3355"].keys():
            references = [PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)]
            disp_references = [PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)]
            for disease in evidence["P3355"][drug].keys():
                if disease != None:
                    disease_qualifier = [PBB_Core.WDItemID(value=disease, prop_nr='P2175', is_qualifier=True)]
                    if "stated_in" in evidence["P3355"][drug][disease].keys():
                        for statedin in evidence["P3355"][drug][disease]["stated_in"]["references"]:
                            references.append(PBB_Core.WDItemID(value=statedin, prop_nr='P248', is_reference=True))
                        for evidence_id in evidence["P3355"][drug][disease]["stated_in"]["id"]:
                            references.append(PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/evidence/"+evidence_id, prop_nr="P854", is_reference=True))
                        prep["P3355"].append(PBB_Core.WDItemID(value=drug, prop_nr='P3355', references=[copy.deepcopy(references)] , qualifiers=copy.deepcopy(disease_qualifier)))

                    if "disputed_by" in evidence["P3355"][drug][disease].keys():
                        for disputedby in evidence["P3355"][drug][disease]["disputed_by"]["references"]:
                            disp_references.append(PBB_Core.WDItemID(value=disputedby, prop_nr='P248', is_reference=True))
                            disease_qualifier.append(PBB_Core.WDItemID(value=disputedby, prop_nr='P1310', is_qualifier=True))
                        for evidence_id in evidence["P3355"][drug][disease]["disputed_by"]["id"]:
                            disp_references.append(PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/evidence/"+evidence_id, prop_nr="P854", is_reference=True))
                        prep["P3355"].append(PBB_Core.WDItemID(value=drug, prop_nr='P3355', references=[copy.deepcopy(disp_references)] , qualifiers=copy.deepcopy(disease_qualifier)))


    if len(evidence["P3356"]) > 0:
        prep["P3356"] = []
        for disease in evidence["P3356"].keys():
            if disease != None:
                references = [PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)]
                disp_references = [PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)]
                if "stated_in" in evidence["P3356"][disease]["references"].keys():
                    references.append(PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/evidence/"+evidence["P3356"][disease]["id"], prop_nr="P854", is_reference=True))
                    for statedin in evidence["P3356"][disease]["references"]["stated_in"]:
                        references.append(PBB_Core.WDItemID(value=statedin, prop_nr='P248', is_reference=True))
                    prep["P3356"].append(PBB_Core.WDItemID(value=disease, prop_nr='P3356', references=[copy.deepcopy(references)]))

                if "disputed_by" in evidence["P3356"][disease]["references"].keys():
                    disp_references.append(PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/evidence/"+evidence["P3356"][disease]["id"], prop_nr="P854", is_reference=True))
                    disputed_qualifier = []
                    for disputedby in evidence["P3356"][disease]["references"]["disputed_by"]:
                        disputed_qualifier.append(PBB_Core.WDItemID(value=disputedby, prop_nr='P1310', is_qualifier=True))
                        disp_references.append(PBB_Core.WDItemID(value=disputedby, prop_nr='P248', is_reference=True))
                    prep["P3356"].append(PBB_Core.WDItemID(value=disease, prop_nr='P3356', references=[copy.deepcopy(disp_references)], qualifiers=copy.deepcopy(disputed_qualifier)))


    if len(evidence["P3357"]) > 0:
        prep["P3357"] = []

        for disease in evidence["P3357"].keys():
            if disease != None:
                references = [PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)]
                disp_references = [PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)]

                if "stated_in" in evidence["P3357"][disease]["references"].keys():
                    references.append(PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/evidence/"+evidence["P3357"][disease]["id"], prop_nr="P854", is_reference=True))
                    for statedin in evidence["P3357"][disease]["references"]["stated_in"]:
                        references.append(PBB_Core.WDItemID(value=statedin, prop_nr='P248', is_reference=True))
                    prep["P3357"].append(PBB_Core.WDItemID(value=disease, prop_nr='P3357', references=[copy.deepcopy(references)]))

                if "disputed_by" in evidence["P3357"][disease]["references"].keys():
                    disputed_qualifier = []
                    disp_references.append(PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/evidence/"+evidence["P3357"][disease]["id"], prop_nr="P854", is_reference=True))
                    for disputedby in evidence["P3357"][disease]["references"]["disputed_by"]:
                        disp_references.append(PBB_Core.WDItemID(value=disputedby, prop_nr='P248', is_reference=True))
                        disputed_qualifier.append(PBB_Core.WDItemID(value=disputedby, prop_nr='P1310', is_qualifier=True))
                    prep["P3357"].append(PBB_Core.WDItemID(value=disease, prop_nr='P3357', references=[copy.deepcopy(disp_references)], qualifiers=copy.deepcopy(disputed_qualifier)))

    if len(evidence["P3358"]) > 0:
        prep["P3358"] = []

        for disease in evidence["P3358"].keys():
            if disease != None:
                references = [PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)]
                disp_references = [PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)]

                if "stated_in" in evidence["P3358"][disease]["references"].keys():
                    references.append(PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/evidence/"+evidence["P3358"][disease]["id"], prop_nr="P854", is_reference=True))
                    for statedin in evidence["P3358"][disease]["references"]["stated_in"]:
                        references.append(PBB_Core.WDItemID(value=statedin, prop_nr='P248', is_reference=True))
                    prep["P3358"].append(PBB_Core.WDItemID(value=disease, prop_nr='P3358', references=[copy.deepcopy(references)]))
                if "disputed_by" in evidence["P3358"][disease]["references"].keys():
                    disputed_qualifier = []
                    disp_references.append(PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/evidence/"+evidence["P3358"][disease]["id"], prop_nr="P854", is_reference=True))
                    for disputedby in evidence["P3358"][disease]["references"]["disputed_by"]:
                        disputed_qualifier.append(PBB_Core.WDItemID(value=disputedby, prop_nr='P1310', is_qualifier=True))
                        disp_references.append(PBB_Core.WDItemID(value=disputedby, prop_nr='P248', is_reference=True))
                    prep["P3358"].append(PBB_Core.WDItemID(value=disease, prop_nr='P3358', references=[copy.deepcopy(disp_references)], qualifiers=copy.deepcopy(disputed_qualifier)))

    if len(evidence["P3359"]) > 0:
        prep["P3359"] = []
        timeStringNow = strftime("+%Y-%m-%dT00:00:00Z", gmtime())

        for disease in evidence["P3359"].keys():
            if disease != None:
                references = [PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)]
                disp_references = [PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)]

                if "stated_in" in evidence["P3359"][disease]["references"].keys():
                    references.append(PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/evidence/"+evidence["P3359"][disease]["id"], prop_nr="P854", is_reference=True))
                    for statedin in evidence["P3359"][disease]["references"]["stated_in"]:
                        references.append(PBB_Core.WDItemID(value=statedin, prop_nr='P248', is_reference=True))
                    prep["P3359"].append(PBB_Core.WDItemID(value=disease, prop_nr='P3359', references=[copy.deepcopy(references)]))
                if "disputed_by" in evidence["P3359"][disease]["references"].keys():
                    disputed_qualifier = []
                    disp_references.append(PBB_Core.WDUrl("https://civic.genome.wustl.edu/links/evidence/"+evidence["P3359"][disease]["id"], prop_nr="P854", is_reference=True))
                    for disputedby in evidence["P3359"][disease]["references"]["disputed_by"]:
                        disputed_qualifier.append(PBB_Core.WDItemID(value=disputedby, prop_nr='P1310', is_qualifier=True))
                        disp_references.append(PBB_Core.WDItemID(value=disputedby, prop_nr='P248', is_reference=True))
                    prep["P3359"].append(PBB_Core.WDItemID(value=disease, prop_nr='P3359', references=[copy.deepcopy(disp_references)], qualifiers=copy.deepcopy(disputed_qualifier)))

    data2add = []
    for key in prep.keys():
        for statement in prep[key]:
            data2add.append(statement)
            print(statement.prop_nr, statement.value)

    pprint(prep)
    name = variant_data["name"]
    wdPage = PBB_Core.WDItemEngine( item_name=name, data=data2add, server="www.wikidata.org", domain="genes", fast_run=fast_run, fast_run_base_filter=fast_run_base_filter)
    synonyms = []
    if name not in ignore_synonym_list:
        synonyms.append(name)
    else:
        if name == "EXPRESSION":
            wdPage.set_label("expressie van "+variant_data["entrez_name"], "nl")
        elif name == "BIALLELIC INACTIVATION":
            wdPage.set_label("biallelische inactivatie van "+variant_data["entrez_name"], "nl")



    wdPage.set_label(variant_data["entrez_name"]+" "+ name, "en")
    # wdPage.set_label(variant_data["entrez_name"]+" "+ name, "nl")


    if wdPage.get_description(lang='en') == "":
        wdPage.set_description("genetic variant", "en")
    if wdPage.get_description(lang='nl') == "":
        wdPage.set_description("gen variant", "nl")
    if len(variant_data["variant_aliases"])>0:
        for alias in variant_data["variant_aliases"]:
            synonyms.append(alias)
    if len(synonyms) > 0:
        wdPage.set_aliases(aliases=synonyms, lang='en', append=True)

    wd_json_representation = wdPage.get_wd_json_representation()

    pprint(wd_json_representation)

    print(wdPage.write(logincreds))

 # except Exception as e:
    """
                print(traceback.format_exc())
                PBB_Core.WDItemEngine.log('ERROR', '{main_data_id}, "{exception_type}", "{message}", {wd_id}, {duration}'.format(
                        main_data_id=variant_id,
                        exception_type=type(e),
                        message=e.__str__(),
                        wd_id='-',
                        duration=time.time()
                    ))"" \
    """
