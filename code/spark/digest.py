
# coding: utf-8

# In[15]:

from pyspark import SparkContext, SparkConf
from pyspark.sql import SQLContext
from pyspark.sql import HiveContext
import json
import time
import sys
from  scipy.stats import fisher_exact, ttest_ind

content = [line.rstrip() for line in open('jobsArguments.conf')]

analysisName=content[0]
scope=content[1]
scale=content[2]
sqlControl=content[3]
sqlCase=content[4]
group1name=content[5]
group2name=content[6]
controlMAF=content[7]
caseMAF=content[8]
pathVariants=content[9]

#nPartitions=8
conf = (SparkConf()
#         .setMaster("local["+str(nPartitions)+"]")
         .setAppName(analysisName)
       )
#sc.stop()
sc = SparkContext(conf=conf)


#parquetFile = sqlContext.read.parquet("/user/hive/warehouse/gvr4.db/variantsulb")
#parquetFile = sqlContext.read.parquet("/Users/yalb/Projects/Github/Docker/cdh54_4_add1000g/variants2")
#parquetFile = sqlContext.read.parquet("hdfs://127.0.0.1:8020/user/hive/warehouse/gvr.db/test")
#parquetFile = sqlContext.read.parquet("hdfs://localhost/user/hive/warehouse/gvr3.db/variants")



# In[16]:

#sqlContext = HiveContext(sc) #sqlContext._get_hive_ctx() #HiveContext(sc) 
sqlContext = SQLContext(sc)
sqlContext.sql("SET spark.sql.parquet.binaryAsString=true")

parquetFile = sqlContext.read.parquet(pathVariants)
parquetFile.registerTempTable("variantData");


# In[17]:

#Input is vector patient, chr, pos, ref, alt, gene_symbol, zygosity
def createKey_VariantGene(variantData):
    #ID is chr:pos:ref:alt
    ID=variantData[1]+":"+str(variantData[2])+":"+variantData[3]+":"+variantData[4]
    
    #return ID, gene_symbol, patient, zygosity
    zygosity=1
    if variantData[6]=="Homozygous":
    #if variantData[6]==2:
        zygosity=2
    patientsID_dictionnary=patientsID_dictionnary_b.value
    patientIndex=patientsID_dictionnary[variantData[0]]
    return ((ID,variantData[5]),(patientIndex,zygosity))

#variantGeneEntry: key is (variantID,gene), value is (patientIndex,zygosity)
def geneAsKey(variantGeneEntry):    
    return (variantGeneEntry[0][1],(variantGeneEntry[0][0],variantGeneEntry[1]))

def createPairs(k,v,idList):
    idListOthers=idList[:]
    idListOthers.remove(k)
    return [(str(sorted([k,idElt])),(sorted([k,idElt]),v)) for idElt in idListOthers]

def getVariantID(key_VariantGene):
    return key_VariantGene[0]


# In[18]:

def getGenotypeVector(genotypeList):
    genotypeVector=[0]*len(patientsID_dictionnary_b.value)
    if len(genotypeList)>0:
        for j in range(0,len(genotypeList)):
            genotypeVector[genotypeList[j][0]]=genotypeList[j][1]
        
        sumCase=float(sum([int(x>0) for x in genotypeVector[0:patientsID_split_index_b.value]]))
        sumControl=float(sum([int(x>0) for x in genotypeVector[patientsID_split_index_b.value:len(patientsID_dictionnary_b.value)]]))
    
        ratioCase=sumCase/patientsID_split_index_b.value
        ratioControl=sumControl/(len(patientsID_dictionnary_b.value)-patientsID_split_index_b.value)
        
        if (ratioCase>float(caseMAF_b.value)) or (ratioControl>float(controlMAF_b.value)):
            genotypeVector=[0]*len(patientsID_dictionnary_b.value)
        
    return genotypeVector        
    


# In[67]:

def getGenotypeVectorByGene(geneID,variantList):
    genotypeVectorByGene=[0]*len(patientsID_dictionnary_b.value)
    
    if len(variantList)>0:
        #Go through list of variants
        for i in range(0,len(variantList)):
            #Get variant ID, and list of sample_index,genotype
            (variantID,genotypeList)=variantList[i]
            if genotypeList.__class__==tuple:
                genotypeList=[genotypeList]
            else:
                genotypeList=list(genotypeList)
            
            #Get genotype vector for current variantID
            genotypeVector=getGenotypeVector(genotypeList)
            #And sum with previous genotype vectors
            genotypeVectorByGene=[x+y for x,y in zip(genotypeVectorByGene,genotypeVector)]
    
    return genotypeVectorByGene


# In[20]:

#variantList is [(locusID,[genotype])]
def scoreVariant(key_VariantGene,value_GenotypeList):
    genotypeList=list(value_GenotypeList)
    
    patientsID_dictionnary=patientsID_dictionnary_b.value
    patientsID_split_index=patientsID_split_index_b.value
    
    genotypeVector=getGenotypeVector(genotypeList)
    
    #sumCase=float(sum([int(x>0) for x in genotypeVector[0:patientsID_split_index]]))
    #sumControl=float(sum([int(x>0) for x in genotypeVector[patientsID_split_index:len(patientsID_dictionnary)]]))
    sumCase=float(sum([x for x in genotypeVector[0:patientsID_split_index]]))
    sumControl=float(sum([x for x in genotypeVector[patientsID_split_index:len(patientsID_dictionnary)]]))
    
    ratioCase=sumCase/patientsID_split_index
    ratioControl=sumControl/(len(patientsID_dictionnary)-patientsID_split_index)
        
    score=ratioCase-ratioControl
    #pvalue=fisher_exact([[sumCase,patientsID_split_index-sumCase],[sumControl,len(patientsID_dictionnary)-patientsID_split_index]],'greater')[1]
    pvalue=ttest_ind(genotypeVector[0:patientsID_split_index],genotypeVector[patientsID_split_index:len(patientsID_dictionnary)])[1]/2
    
    
    #if score>0:
    return (key_VariantGene,(score,pvalue,ratioCase,ratioControl,sumCase,sumControl))


# In[21]:

def scoreVariantPair(variantIDpair,value_GenotypeListPair):
    
    genotypeListPair=list(value_GenotypeListPair)
    
    patientsID_dictionnary=patientsID_dictionnary_b.value
    patientsID_split_index=patientsID_split_index_b.value
    
    score=0
    if len(genotypeListPair)==2:
        (variantID,genotypeList1)=genotypeListPair[0]
        (variantID,genotypeList2)=genotypeListPair[1]
        
        variantID1=variantID[0]
        variantID2=variantID[1]
        
        genotypeList1=list(genotypeList1)
        genotypeList2=list(genotypeList2)
        
        genotypeVector1=getGenotypeVector(genotypeList1)
        genotypeVector2=getGenotypeVector(genotypeList2)
        
        genotypeVector=[int(x>0 and y>0) for x,y in zip(genotypeVector1,genotypeVector2)]
        
        sumCase=float(sum([int(x>0) for x in genotypeVector[0:patientsID_split_index]]))
        ratioCase=sumCase/patientsID_split_index
        sumControl=float(sum([int(x>0) for x in genotypeVector[(patientsID_split_index+1):len(patientsID_dictionnary)]]))
        ratioControl=sumControl/(len(patientsID_dictionnary)-patientsID_split_index)
        
        score=ratioCase-ratioControl
        pvalue=fisher_exact([[sumCase,patientsID_split_index-sumCase],[sumControl,len(patientsID_dictionnary)-patientsID_split_index]],'greater')[1]
        
        #if score>0:
        return (variantIDpair,((variantID1,variantID2),score,pvalue,ratioCase,ratioControl,sumCase,sumControl))



# In[22]:

#variantList is [(locusID,[sample_index,genotype])]
def scoreGene(geneID,variantList):
    variantList=list(variantList)
    
    patientsID_dictionnary=patientsID_dictionnary_b.value
    genoSum=[0]*len(patientsID_dictionnary)
    
    patientsID_split_index=patientsID_split_index_b.value
    
    genotypeVectorByGene=getGenotypeVectorByGene(geneID,variantList)
    
    sumCase=float(sum([int(x>0) for x in genotypeVectorByGene[0:patientsID_split_index]]))
    sumControl=float(sum([int(x>0) for x in genotypeVectorByGene[patientsID_split_index:len(patientsID_dictionnary)]]))
    
    ratioCase=sumCase/patientsID_split_index
    ratioControl=sumControl/(len(patientsID_dictionnary)-patientsID_split_index)
        
    score=ratioCase-ratioControl
    pvalue=fisher_exact([[sumCase,patientsID_split_index-sumCase],[sumControl,len(patientsID_dictionnary)-patientsID_split_index]],'greater')[1]
    #pvalue=ttest_ind(genotypeVectorByGene[0:patientsID_split_index],genotypeVectorByGene[patientsID_split_index:len(patientsID_dictionnary)])[1]/2
        
    #if score>0:
    return (geneID,(score,pvalue,ratioCase,ratioControl,sumCase,sumControl))


# In[23]:

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
        
        genoSum1=getGenotypeVectorByGene(gene1,variantList1)
        genoSum2=getGenotypeVectorByGene(gene2,variantList2)
        
        genoSum=[int(x>0 and y>0) for x,y in zip(genoSum1,genoSum2)]
        sumCase=float(sum([int(x>0) for x in genoSum[0:patientsID_split_index]]))
        sumControl=float(sum([int(x>0) for x in genoSum[(patientsID_split_index+1):len(patientsID_dictionnary)]]))
        
        ratioCase=sumCase/patientsID_split_index
        ratioControl=sumControl/(len(patientsID_dictionnary)-patientsID_split_index)
        
        score=ratioCase-ratioControl
        pvalue=fisher_exact([[sumCase,patientsID_split_index-sumCase],[sumControl,len(patientsID_dictionnary)-patientsID_split_index]],'greater')[1]
        
        #if score>0:
        return (gene_symbol_pair,((gene1,gene2),score,pvalue,ratioCase,ratioControl,sumCase,sumControl))



# In[24]:

start_time = time.time()

variants_case = sqlContext.sql("SELECT sample_id,chr,pos,ref,alt,gene_symbol,zygosity FROM variantData "+sqlCase)
patientsID_case=sorted(variants_case.map(lambda v:v[0]).distinct().collect())

variants_control= sqlContext.sql("SELECT sample_id,chr,pos,ref,alt,gene_symbol,zygosity FROM variantData "+sqlControl)
patientsID_control=sorted(variants_control.map(lambda v:v[0]).distinct().collect())

patientsID=patientsID_case+patientsID_control
patientsID_dictionnary=dict(zip(patientsID,range(len(patientsID))))

patientsID_split_index_b = sc.broadcast(len(patientsID_case))
patientsID_dictionnary_b = sc.broadcast(patientsID_dictionnary)

variants=variants_control.unionAll(variants_case)
variants_grouped=variants.map(createKey_VariantGene)

controlMAF_b=sc.broadcast(controlMAF)
caseMAF_b=sc.broadcast(caseMAF)

#variants_grouped.count()


# In[68]:

#start_time = time.time()
ntests=0

if scope=='monogenic':
    if scale=='variant':
        ntests=variants_grouped.count()
        scores=variants_grouped.map(lambda (k,v):scoreVariant(k,v)).filter(lambda x:x is not None).takeOrdered(1000, key=lambda (k,(v1,v2,v3,v4,v5,v6)): -v1)

    if scale=='gene':
        variants_grouped_by_gene=variants_grouped.map(geneAsKey).groupByKey()
        ntests=variants_grouped_by_gene.count()
        scores=variants_grouped_by_gene.map(lambda (k,v):scoreGene(k,v)).filter(lambda x:x is not None).takeOrdered(1000, key=lambda (k,(v1,v2,v3,v4,v5,v6)): -v1)

if scope=='digenic':
    if scale=='variant':
        variantsID=variants_grouped.keys().map(getVariantID).collect()
        variants_grouped_by_pairs=variants_grouped.flatMap(lambda (k,v):createPairs(k[0],v,variantsID)).groupByKey()
        scores=variants_grouped_by_pairs.map(lambda (k,v):scoreVariantPair(k,v)).filter(lambda x:x is not None).takeOrdered(1000, key=lambda (k,(variants,v1,v2,v3,v4,v5,v6)): -v1)
        ntests=len(variantsID)*(len(variantsID)+1)/2
   
    if scale=='gene':
        variants_grouped_by_gene=variants_grouped.map(geneAsKey).groupByKey()
        genesID=variants_grouped_by_gene.keys().collect()
        variants_grouped_by_gene_pairs=variants_grouped_by_gene.flatMap(lambda (k,v):createPairs(k,v,genesID)).groupByKey()
        scores=variants_grouped_by_gene_pairs.map(lambda (k,v):scoreGenePair(k,v)).filter(lambda x:x is not None).takeOrdered(1000, key=lambda (k,(genes,v1,v2,v3,v4,v5,v6)): -v1)
        ntests=len(genesID)*(len(genesID)+1)/2
    
end_time=time.time()
runtime=end_time - start_time
print(runtime)


# In[69]:

len(scores)


# In[70]:

scores[0:3]


# In[14]:

scores=[analysisName,scale,scope,start_time,end_time,runtime,scores,patientsID_case,patientsID_control,group1name,group2name,ntests]

with open(analysisName+'.txt', 'w') as outfile:
    json.dump(scores, outfile)
    


# In[ ]:

sc.stop()

