
# coding: utf-8

# In[88]:

from pyspark import SparkContext, SparkConf
from pyspark.sql import SQLContext
from pyspark.sql import HiveContext
import json
import time
import sys
from  scipy.stats import fisher_exact

content = [line.rstrip() for line in open('jobsArguments.conf')]

analysisName=content[0]
scope=content[1]
scale=content[2]
sqlControl=content[3]
sqlCase=content[4]
group1name=content[5]
group2name=content[6]

nPartitions=16
conf = (SparkConf()
         .setMaster("local["+str(nPartitions)+"]")
         .setAppName(analysisName)
#         .set("spark.executor.memory", "5g")
#         .set("spark.driver.memory", "5g")
#         .set("spark.python.worker.memory", "5g")
       )
#sc.stop()
sc = SparkContext(conf=conf)


#parquetFile = sqlContext.read.parquet("/user/hive/warehouse/gvr4.db/variantsulb")
#parquetFile = sqlContext.read.parquet("/Users/yalb/Projects/Github/Docker/cdh54_4_add1000g/variants2")
#parquetFile = sqlContext.read.parquet("hdfs://127.0.0.1:8020/user/hive/warehouse/gvr.db/test")
#parquetFile = sqlContext.read.parquet("hdfs://localhost/user/hive/warehouse/gvr3.db/variants")



# In[98]:

#sqlContext = HiveContext(sc) #sqlContext._get_hive_ctx() #HiveContext(sc) 
sqlContext = SQLContext(sc)
sqlContext.sql("SET spark.sql.parquet.binaryAsString=true")

#parquetFile = sqlContext.read.parquet("hdfs://node001:8020/user/hive/warehouse/highlander.db/exomes_hc")

parquetFile = sqlContext.read.parquet("/Users/yalb/Projects/Github/digest/variantsulb")
parquetFile.registerTempTable("parquetFile");


# In[140]:

#analysisName="control_vs_neurodev_rare_digenic"
#group1name="control_ulb_rare_damaging"
#group2name="neurodev_ulb_rare_damaging"
scope="monogenic"
#scale="gene"
#controlMAF=0.5


# In[100]:

#RDDtest = sqlContext.sql("SELECT distinct patient from parquetFile")


# In[101]:

#RDDtest.count()


# In[102]:

#Input is vector patient, chr, pos, ref, alt, gene_symbol, zygosity
def createKey_VariantGene(variantData):
    #ID is chr:pos:ref:alt
    ID=variantData[1]+":"+str(variantData[2])+":"+variantData[3]+":"+variantData[4]
    #return ID, gene_symbol, patient, zygosity
    zygosity=2
    if variantData[6]=="Heterozygous":
        zygosity=1
    patientsID_dictionnary=patientsID_dictionnary_b.value
    patientIndex=patientsID_dictionnary[variantData[0]]
    return ((ID,variantData[5]),(patientIndex,zygosity))

def buildVariantVector(ID,variantData,patientsID):
    variantData=list(variantData)
    genotype=[0]*len(patientsID)
    
    #Get sampleID/Genotype for each variant
    for i in range(0,len(variantData)):
        genotype[variantData[i][1]]=variantData[i][2]
    
    return ((ID,variantData[0][0]),genotype)



# In[103]:

#variantGeneEntry: key is (variantID,gene), value is (patientIndex,zygosity)
def geneAsKey(variantGeneEntry):    
    return (variantGeneEntry[0][1],(variantGeneEntry[0][0],variantGeneEntry[1]))

def makePairParts(k,v,nbPart):
    result=[]
    for i in range(0,nbPart):
        result.append(((k,i),v))
        
    return [(str(sorted([k,i])),(v)) for i in range(0,nbPart)]

def f(splitIndex ,v): 
    return [(splitIndex,list(v))]


# In[104]:

def scoreVariantUnivariate(k,variantData):
    variantData=list(variantData)
    
    score=0
    sumControl=0
    
    sumCase=sum([int(x>0) for x in variantData[0]])
    sumControl=sum([int(x>0) for x in variantData[1]])
    
    score=sumCase#-sumControl
    if sumControl>0:
        score=0
        
    if score>0:
        return (k,(score,sumCase,sumControl))


# In[105]:

def getGenotypeVectorByGene(gene_symbol,variantList,patientsID_dictionnary,patientsID_split_index):
    genoSum=[0]*len(patientsID_dictionnary)
    
    if len(variantList)>0:
        #Go through list of variants
        for i in range(0,len(variantList)):
            #Get variant ID, and list of sample_index,genotype
            (variantID,sample_geno_list)=variantList[i]
            sample_geno_list=list(sample_geno_list)
            
            #Go through list of sample_index,genotype
            for j in range(0,len(sample_geno_list)):
                genoSum[sample_geno_list[j][0]]=genoSum[sample_geno_list[j][0]]+sample_geno_list[j][1]
    return genoSum


# In[146]:

#variantList is [(locusID,[sample_index,genotype])]
def scoreGene(gene_symbol,variantList):
    variantList=list(variantList)
    
    patientsID_dictionnary=patientsID_dictionnary_b.value
    genoSum=[0]*len(patientsID_dictionnary)
    
    patientsID_split_index=patientsID_split_index_b.value
    
    genoSum=getGenotypeVectorByGene(gene_symbol,variantList,patientsID_dictionnary,patientsID_split_index)
    
    sumCase=float(sum([int(x>0) for x in genoSum[0:patientsID_split_index]]))
    sumControl=float(sum([int(x>0) for x in genoSum[(patientsID_split_index+1):len(patientsID_dictionnary)]]))
    
    ratioCase=sumCase/patientsID_split_index
    ratioControl=sumControl/(len(patientsID_dictionnary)-patientsID_split_index)
        
    score=ratioCase-ratioControl
    pvalue=fisher_exact([[sumCase,patientsID_split_index-sumCase],[sumControl,len(patientsID_dictionnary)-patientsID_split_index]],'greater')[1]
        
    if score>0:
        return (gene_symbol,(score,pvalue,ratioCase,ratioControl,sumCase,sumControl))


# In[129]:

def scoreGenePair(gene_symbol_pair,variantList):
    
    variantList=list(variantList)
    
    patientsID_dictionnary=patientsID_dictionnary_b.value
    
    patientsID_split_index=patientsID_split_index_b.value
    
    score=0
    if len(variantList)==2:
        (genes,variantList1)=variantList[0]
        (genes,variantList2)=variantList[1]
        
        gene1=genes[0]
        gene2=genes[1]
        
        variantList1=list(variantList1)
        variantList2=list(variantList2)
        
        genoSum1=getGenotypeVectorByGene(gene1,variantList1,patientsID_dictionnary,patientsID_split_index)
        genoSum2=getGenotypeVectorByGene(gene2,variantList2,patientsID_dictionnary,patientsID_split_index)
        
        genoSum=[int(x>0 and y>0) for x,y in zip(genoSum1,genoSum2)]
        
        sumCase=float(sum([int(x>0) for x in genoSum[0:patientsID_split_index]]))
        ratioCase=sumCase/patientsID_split_index
        sumControl=float(sum([int(x>0) for x in genoSum[(patientsID_split_index+1):len(patientsID_dictionnary)]]))
        ratioControl=sumControl/(len(patientsID_dictionnary)-patientsID_split_index)
        
        score=ratioCase-ratioControl
        pvalue=fisher_exact([[sumCase,patientsID_split_index-sumCase],[sumControl,len(patientsID_dictionnary)-patientsID_split_index]],'greater')[1]
        
        if score>0:
            return (gene_symbol_pair,((gene1,gene2),score,pvalue,ratioCase,ratioControl,sumCase,sumControl))

#Key is (variantID, gene)
def getGene(variantGene_key):
    gene=variantGene_key[1]
    
    return (gene)

def createPairsGenes(k,v,genes):
    return [(str(sorted([k,gene])),(sorted([k,gene]),v)) for gene in genes]

def fillMissing(k,v):
    v=list(v)
    if v[0] is None:
        v[0]=[0]*len(sample_id_case_b.value)
    if v[1] is None:
        v[1]=[0]*len(sample_id_control_b.value)
        
    return (k,v)

def fillMissing(k,v):
    v=list(v)
    if v[1] is None:
        v[1]=[0]*len(dict_patient_control_b.value)
        
    return (k,v)


# In[108]:

float(15)/2


# In[141]:

start_time = time.time()

variants_case = sqlContext.sql("SELECT patient,chr,pos,reference,alternative,gene_symbol,zygosity FROM parquetFile "+sqlCase)
patientsID_case=sorted(variants_case.map(lambda v:v[0]).distinct().collect())

if sqlControl!="NULL":
    variants_control= sqlContext.sql("SELECT patient,chr,pos,reference,alternative,gene_symbol,zygosity FROM parquetFile "+sqlControl)
#    controlMAF=float(controlMAF)
else:
    variants_control=sc.emptyRDD()
#    controlMAF=0   
patientsID_control=sorted(variants_control.map(lambda v:v[0]).distinct().collect())

patientsID=patientsID_case+patientsID_control
patientsID_dictionnary=dict(zip(patientsID,range(len(patientsID))))

patientsID_split_index_b=sc.broadcast(len(patientsID_case))

patientsID_dictionnary_b = sc.broadcast(patientsID_dictionnary)

variants=variants_control.unionAll(variants_case)

variants_grouped=variants.map(createKey_VariantGene).groupByKey()


# In[150]:

#start_time = time.time()
ntests=0

if scope=='monogenic':
    if scale=='variant':
        scores=genoMat.map(lambda (k,v):scoreVariantUnivariate(k,v)).filter(lambda x:x is not None).takeOrdered(10000000, key=lambda (k,(v1,v2,v3)): -v1)

    if scale=='gene':
        variants_grouped_by_gene=variants_grouped.map(geneAsKey).groupByKey()
        ntests=variants_grouped_by_gene.count()
        scores=variants_grouped_by_gene.map(lambda (k,v):scoreGene(k,v)).filter(lambda x:x is not None).takeOrdered(1000, key=lambda (k,(v1,v2,v3,v4,v5,v6)): -v1)
if scope=='digenic':
    variants_grouped_by_gene=variants_grouped.map(geneAsKey).groupByKey()#.flatMap(lambda (k,v):createPairsGenes(k,v,genes)).groupByKey().map(lambda (k,v):scoreGenePair(k,v))#.filter(lambda x:x is not None).takeOrdered(1000, key=lambda (k,(genes,v1,v2,v3)): -v1)
    genes=variants_grouped_by_gene.keys().collect()
    scores=variants_grouped_by_gene.flatMap(lambda (k,v):createPairsGenes(k,v,genes)).groupByKey().map(lambda (k,v):scoreGenePair(k,v)).filter(lambda x:x is not None).takeOrdered(1000, key=lambda (k,(genes,v1,v2,v3,v4,v5,v6)): -v1)
    ntests=len(genes)*(len(genes)+1)/2
    
end_time=time.time()
runtime=end_time - start_time
print(runtime)


# In[151]:

len(scores)


# In[149]:

scores[0:3]


# In[50]:

#analysisName="neurodev_vs_control_digenic_rare_high"


# In[152]:

scores=[analysisName,scale,scope,start_time,end_time,runtime,scores,patientsID_case,patientsID_control,group1name,group2name,ntests]

with open(analysisName+'.txt', 'w') as outfile:
    json.dump(scores, outfile)
    


# In[ ]:

sc.stop()

