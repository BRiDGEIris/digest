
# coding: utf-8

# In[21]:

from pyspark import SparkContext, SparkConf
from pyspark.sql import SQLContext
from pyspark.sql import HiveContext
import json
import time
import sys

content = [line.rstrip() for line in open('jobsArguments.conf')]

analysisName=content[0]
scope=content[1]
scale=content[2]
sqlControl=content[3]
sqlCase=content[4]
group1name=content[5]
group2name=content[6]

nPartitions=4
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



# In[22]:

#sqlContext = HiveContext(sc) #sqlContext._get_hive_ctx() #HiveContext(sc) 
sqlContext = SQLContext(sc)
sqlContext.sql("SET spark.sql.parquet.binaryAsString=true")

#parquetFile = sqlContext.read.parquet("hdfs://node001:8020/user/hive/warehouse/highlander.db/exomes_hc")

parquetFile = sqlContext.read.parquet("/Users/yalb/Projects/Github/digest/variantsulb")
parquetFile.registerTempTable("parquetFile");


# In[3]:

#analysisName="control_vs_neurodev_rare_digenic"
#group1name="control_ulb_rare_damaging"
#group2name="neurodev_ulb_rare_damaging"
#scope="digenic"
#scale="gene"
#controlMAF=0.5


# In[4]:

#RDDtest = sqlContext.sql("SELECT distinct patient from parquetFile")


# In[5]:

#RDDtest.count()


# In[23]:

#Input is vector patient, chr, pos, ref, alt, gene_symbol, zygosity
def createKey_VariantGene(variantData):
    #ID is chr:pos:ref:alt
    ID=variantData[1]+":"+str(variantData[2])+":"+variantData[3]+":"+variantData[4]
    #return ID, gene_symbol, patient, zygosity
    zygosity=2
    if variantData[6]==1:
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



# In[24]:

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


# In[25]:

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


# In[26]:

#variantList is [(locusID,[sample_index,genotype])]
def scoreGeneUnivariate(gene_symbol,variantList):
    variantList=list(variantList)
    
    patientsID_dictionnary=patientsID_dictionnary_b.value
    genoSum=[0]*len(patientsID_dictionnary)
    
    patientsID_split_index=patientsID_split_index_b.value
    
    if len(variantList)>0:
        #Go through list of variants
        for i in range(0,len(variantList)):
            #Get variant ID, and list of sample_index,genotype
            (variantID,sample_geno_list)=variantList[i]
            sample_geno_list=list(sample_geno_list)
            
            #Go through list of sample_index,genotype
            for j in range(0,len(sample_geno_list)):
                genoSum[sample_geno_list[j][0]]=genoSum[sample_geno_list[j][0]]+sample_geno_list[j][1]
        sumCase=sum([int(x>0) for x in genoSum[0:patientsID_split_index]])
        sumControl=sum([int(x>0) for x in genoSum[(patientsID_split_index+1):len(patientsID_dictionnary)]])
        score=sumCase-sumControl

    if score>0:
        return (gene_symbol,(score,sumCase,sumControl))
    #return (gene_symbol,(genoSum))


# In[27]:

def scoreDigenicGene(variantLists):
    variantLists=list(variantLists)
    geno1sumcase=[]
    geno1sumcontrol=[]
    geno2sumcase=[]
    geno2sumcontrol=[]
    score=0
    if len(variantLists)==2:
        (genes,variantList1)=list(variantLists[0])
        (genes,variantList2)=list(variantLists[1])
        gene1=genes[0]
        gene2=genes[1]
        variantList1=list(variantList1)
        variantList2=list(variantList2)
        for i in range(0,len(variantList1)):
            (locus1,geno1)=variantList1[i]
            if geno1sumcase==[]:
                geno1sumcase=[int(x) for x in geno1[0]]
                geno1sumcontrol=[int(x) for x in geno1[1]]
            else:
                geno1sumcase=[int(x)+int(y) for x,y in zip(geno1sumcase,geno1[0])]
                geno1sumcontrol=[int(x)+int(y) for x,y in zip(geno1sumcontrol,geno1[1])]
                
        for i in range(0,len(variantList2)):
            (locus2,geno2)=variantList2[i]
            if geno2sumcase==[]:
                geno2sumcase=[int(x) for x in geno2[0]]
                geno2sumcontrol=[int(x) for x in geno2[1]]
            else:
                geno2sumcase=[int(x)+int(y) for x,y in zip(geno2sumcase,geno2[0])]
                geno2sumcontrol=[int(x)+int(y) for x,y in zip(geno2sumcontrol,geno2[1])]
                
        genosumcase=[int((x>0) & (y>0)) for x,y in zip(geno1sumcase,geno2sumcase)]
        genosumcontrol=[int((x>0) & (y>0)) for x,y in zip(geno1sumcontrol,geno2sumcontrol)]
        
        sumCase=sum([int(x>0) for x in genosumcase])
        sumControl=sum([int(x>0) for x in genosumcontrol])
        
        meanControl=0
        if len(genosumcontrol)>0:
            meanControl=round(float(sumControl)/len(genosumcontrol),2)
        
        score=sumCase
        
        if score>0:
            if meanControl<=controlMAF_b.value:
                return (k,((gene1,gene2),score,sumCase,meanControl))

def getGene(variantData):
    variantGene=variantData[0][1]
    
    return (variantGene)

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


# In[28]:

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


# In[29]:

#start_time = time.time()

if scope=='monogenic':
    if scale=='variant':
        scores=genoMat.map(lambda (k,v):scoreVariantUnivariate(k,v)).filter(lambda x:x is not None).takeOrdered(10000000, key=lambda (k,(v1,v2,v3)): -v1)

    if scale=='gene':
        scores=variants_grouped.map(geneAsKey).groupByKey().map(lambda (k,v):scoreGeneUnivariate(k,v)).filter(lambda x:x is not None).takeOrdered(1000, key=lambda (k,(v1,v2,v3)): -v1)
    
if scope=='digenic':
    genes=genoMat.map(getGene).distinct().takeOrdered(100000).flatMap(lambda (k,v):scoreCompound(k,v)).takeOrdered(100000, key=lambda (k,v1,v2,v3): -v1)
    scores=genoMat.map(splitValues).groupByKey().flatMap(lambda (k,v):createPairsGenes(k,v,genes)).groupByKey().map(lambda (k,v):scoreDigenicGene(k,v)).filter(lambda x:x is not None).takeOrdered(1000, key=lambda (k,(genes,v1,v2,v3)): -v1)

end_time=time.time()
runtime=end_time - start_time
print(runtime)


# In[30]:

scores[0:20]


# In[14]:

scores=[analysisName,scale,scope,start_time,end_time,runtime,scores,patientsID_case,patientsID_control,group1name,group2name]

with open(analysisName+'.txt', 'w') as outfile:
    json.dump(scores, outfile)
    


# In[ ]:

sc.stop()

