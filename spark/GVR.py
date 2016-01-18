
# coding: utf-8

# In[177]:

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



# In[178]:

#sqlContext = HiveContext(sc) #sqlContext._get_hive_ctx() #HiveContext(sc) 
sqlContext = SQLContext(sc)
sqlContext.sql("SET spark.sql.parquet.binaryAsString=true")

#parquetFile = sqlContext.read.parquet("/user/hive/warehouse/gvr4.db/variantsulb")
#parquetFile = sqlContext.read.parquet("/Users/yalb/Projects/Github/Docker/cdh54_4_add1000g/variants2")
#parquetFile = sqlContext.read.parquet("hdfs://127.0.0.1:8020/user/hive/warehouse/gvr.db/test")
#parquetFile = sqlContext.read.parquet("hdfs://localhost/user/hive/warehouse/gvr3.db/variants")

parquetFile = sqlContext.read.parquet("/Users/yalb/Projects/Github/digest/variantsulb")
parquetFile.registerTempTable("parquetFile");


# In[179]:

#analysisName="control_vs_neurodev_rare_digenic"
#group1name="control_ulb_rare_damaging"
#group2name="neurodev_ulb_rare_damaging"
#scope="digenic"
#scale="gene"
#controlMAF=0.5


# In[180]:

#RDDtest = sqlContext.sql("SELECT distinct patient from parquetFile")


# In[181]:

#RDDtest.count()


# In[182]:

#Input is vector patient, chr, pos, ref, alt, gene_symbol, zygosity
def splitByVariantID(variantData,patientsID):
    #ID is chr:pos:ref:alt
    ID=variantData[1]+":"+str(variantData[2])+":"+variantData[3]+":"+variantData[4]
    #return ID, gene_symbol, patient, zygosity
    zygosity=2
    if variantData[6]==1:
        zygosity=1
    patientIndex=patientsID[variantData[0]]
    return (ID,(variantData[5],patientIndex,zygosity))

def buildVariantVector(ID,variantData,patientsID):
    variantData=list(variantData)
    genotype=[0]*len(patientsID)
    
    #Get sampleID/Genotype for each variant
    for i in range(0,len(variantData)):
        genotype[variantData[i][1]]=variantData[i][2]
    
    return ((ID,variantData[0][0]),genotype)



# In[183]:

def geneAsKey(variantData):    
    return (variantData[0][1],(variantData[0][0],variantData[1]))

def makePairParts(k,v,nbPart):
    result=[]
    for i in range(0,nbPart):
        result.append(((k,i),v))
        
    return [(str(sorted([k,i])),(v)) for i in range(0,nbPart)]

def f(splitIndex ,v): 
    return [(splitIndex,list(v))]


# In[184]:

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


# In[185]:

#variantList is [(locusID,[genotypeCase,genotypeControl])]
def scoreGeneUnivariate(k,variantList):
    variantList=list(variantList)
    
    genosumcase=[]
    genosumcontrol=[]
    if len(variantList)>0:
        for i in range(0,len(variantList)):
            (locus,geno)=variantList[i]
            if genosumcase==[]:
                genosumcase=[int(x) for x in geno[0]]
                genosumcontrol=[int(x) for x in geno[1]]
            else:
                genosumcase=[int(x)+int(y) for x,y in zip(genosumcase,geno[0])]
                genosumcontrol=[int(x)+int(y) for x,y in zip(genosumcontrol,geno[1])]
                
        sumCase=sum([int(x>0) for x in genosumcase])
        sumControl=sum([int(x>0) for x in genosumcontrol])
        score=sumCase-sumControl

    if score>0:
        return (k,(score,sumCase,sumControl))


# In[186]:

def scoreDigenicGene(k,variantLists):
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


# In[187]:

start_time = time.time()
#sqlCase="sample_id IN('HG00096','HG00097','HG00099','HG00100','HG00101','HG00102','HG00103','HG00105','HG00106','HG00107','HG00108','HG00109','HG00110','HG00111','HG00112','HG00113','HG00114','HG00115','HG00116','HG00117','HG00118','HG00119','HG00120','HG00121','HG00122','HG00123','HG00125','HG00126','HG00127','HG00128','HG00129','HG00130','HG00131','HG00132','HG00133','HG00136','HG00137','HG00138','HG00139','HG00140','HG00141','HG00142','HG00143','HG00145','HG00146','HG00148','HG00149','HG00150','HG00151','HG00154','HG00155','HG00157','HG00158','HG00159','HG00160','HG00171','HG00173','HG00174','HG00176','HG00177','HG00178','HG00179','HG00180','HG00181','HG00182','HG00183','HG00185','HG00186','HG00187','HG00188','HG00189','HG00190','HG00231','HG00232','HG00233','HG00234','HG00235','HG00236','HG00237','HG00238','HG00239','HG00240','HG00242','HG00243','HG00244','HG00245','HG00246','HG00250','HG00251','HG00252','HG00253','HG00254','HG00255','HG00256','HG00257','HG00258','HG00259','HG00260','HG00261','HG00262')"
#sqlControl="sample_id IN('HG00096','HG00097','HG00099','HG00100','HG00101','HG00102','HG00103','HG00105','HG00106','HG00107','HG00108','HG00109','HG00110','HG00111','HG00112','HG00113','HG00114','HG00115','HG00116','HG00117','HG00118','HG00119','HG00120','HG00121','HG00122','HG00123','HG00125','HG00126','HG00127','HG00128','HG00129','HG00130','HG00131','HG00132','HG00133','HG00136','HG00137','HG00138','HG00139','HG00140','HG00141','HG00142','HG00143','HG00145','HG00146','HG00148','HG00149','HG00150','HG00151','HG00154','HG00155','HG00157','HG00158','HG00159','HG00160','HG00171','HG00173','HG00174','HG00176','HG00177','HG00178','HG00179','HG00180','HG00181','HG00182','HG00183','HG00185','HG00186','HG00187','HG00188','HG00189','HG00190','HG00231','HG00232','HG00233','HG00234','HG00235','HG00236','HG00237','HG00238','HG00239','HG00240','HG00242','HG00243','HG00244','HG00245','HG00246','HG00250','HG00251','HG00252','HG00253','HG00254','HG00255','HG00256','HG00257','HG00258','HG00259','HG00260','HG00261','HG00262')"

#sample_id_case=['HG00096','HG00097','HG00099','HG00100','HG00101','HG00102','HG00103','HG00105','HG00106','HG00107','HG00108','HG00109','HG00110','HG00111','HG00112','HG00113','HG00114','HG00115','HG00116','HG00117','HG00118','HG00119','HG00120','HG00121','HG00122','HG00123','HG00125','HG00126','HG00127','HG00128','HG00129','HG00130','HG00131','HG00132','HG00133','HG00136','HG00137','HG00138','HG00139','HG00140','HG00141','HG00142','HG00143','HG00145','HG00146','HG00148','HG00149','HG00150','HG00151','HG00154','HG00155','HG00157','HG00158','HG00159','HG00160','HG00171','HG00173','HG00174','HG00176','HG00177','HG00178','HG00179','HG00180','HG00181','HG00182','HG00183','HG00185','HG00186','HG00187','HG00188','HG00189','HG00190','HG00231','HG00232','HG00233','HG00234','HG00235','HG00236','HG00237','HG00238','HG00239','HG00240','HG00242','HG00243','HG00244','HG00245','HG00246','HG00250','HG00251','HG00252','HG00253','HG00254','HG00255','HG00256','HG00257','HG00258','HG00259','HG00260','HG00261','HG00262']
#sample_id_control=['HG00096','HG00097','HG00099','HG00100','HG00101','HG00102','HG00103','HG00105','HG00106','HG00107','HG00108','HG00109','HG00110','HG00111','HG00112','HG00113','HG00114','HG00115','HG00116','HG00117','HG00118','HG00119','HG00120','HG00121','HG00122','HG00123','HG00125','HG00126','HG00127','HG00128','HG00129','HG00130','HG00131','HG00132','HG00133','HG00136','HG00137','HG00138','HG00139','HG00140','HG00141','HG00142','HG00143','HG00145','HG00146','HG00148','HG00149','HG00150','HG00151','HG00154','HG00155','HG00157','HG00158','HG00159','HG00160','HG00171','HG00173','HG00174','HG00176','HG00177','HG00178','HG00179','HG00180','HG00181','HG00182','HG00183','HG00185','HG00186','HG00187','HG00188','HG00189','HG00190','HG00231','HG00232','HG00233','HG00234','HG00235','HG00236','HG00237','HG00238','HG00239','HG00240','HG00242','HG00243','HG00244','HG00245','HG00246','HG00250','HG00251','HG00252','HG00253','HG00254','HG00255','HG00256','HG00257','HG00258','HG00259','HG00260','HG00261','HG00262']

RDDcase = sqlContext.sql("SELECT patient,chr,pos,reference,alternative,gene_symbol,zygosity FROM parquetFile "+sqlCase)
patient_case=sorted(RDDcase.map(lambda v:v[0]).distinct().collect())
dict_patient_case=dict(zip(patient_case,range(len(patient_case))))
dict_patient_case_b = sc.broadcast(dict_patient_case)

if sqlControl!="NULL":
    RDDcontrol= sqlContext.sql("SELECT patient,chr,pos,reference,alternative,gene_symbol,zygosity FROM parquetFile "+sqlControl)
#    controlMAF=float(controlMAF)
else:
    RDDcontrol=sc.emptyRDD()
#    controlMAF=0   
patient_control=sorted(RDDcontrol.map(lambda v:v[0]).distinct().collect())
dict_patient_control=dict(zip(patient_control,range(len(patient_control))))
dict_patient_control_b = sc.broadcast(dict_patient_control)


#controlMAF_b=sc.broadcast(controlMAF)

genoMatCase=RDDcase.map(lambda v:splitByVariantID(v,dict_patient_case)).groupByKey()
genoMatCase=genoMatCase.map(lambda (k,v):buildVariantVector(k,v,dict_patient_case))

genoMatControl=RDDcontrol.map(lambda v:splitByVariantID(v,dict_patient_control)).groupByKey()
genoMatControl=genoMatControl.map(lambda (k,v):buildVariantVector(k,v,dict_patient_control))

genoMat=genoMatCase.leftOuterJoin(genoMatControl).map(lambda (k,v): fillMissing(k,v))
#genoMat=genoMatCase.fullOuterJoin(genoMatControl).map(lambda (k,v): fillMissing(k,v))


# In[188]:

genoMat.count()


# In[189]:

#start_time = time.time()


if scope=='monogenic':
    if scale=='variant':
        scores=genoMat.map(lambda (k,v):scoreVariantUnivariate(k,v)).filter(lambda x:x is not None).takeOrdered(10000000, key=lambda (k,(v1,v2,v3)): -v1)

    if scale=='gene':
        scores=genoMat.map(geneAsKey).groupByKey().map(lambda (k,v):scoreGeneUnivariate(k,v)).filter(lambda x:x is not None).takeOrdered(1000, key=lambda (k,(v1,v2,v3)): -v1)
    
if scope=='digenic':
    genes=genoMat.map(getGene).distinct().takeOrdered(100000).flatMap(lambda (k,v):scoreCompound(k,v)).takeOrdered(100000, key=lambda (k,v1,v2,v3): -v1)
    scores=genoMat.map(splitValues).groupByKey().flatMap(lambda (k,v):createPairsGenes(k,v,genes)).groupByKey().map(lambda (k,v):scoreDigenicGene(k,v)).filter(lambda x:x is not None).takeOrdered(1000, key=lambda (k,(genes,v1,v2,v3)): -v1)

end_time=time.time()
runtime=end_time - start_time
print(runtime)


# In[190]:

#scores


# In[192]:

scores=[analysisName,scale,scope,start_time,end_time,runtime,scores,patient_case,patient_control,group1name,group2name]

with open("analyses/"+analysisName+'.txt', 'w') as outfile:
    json.dump(scores, outfile)
    


# In[ ]:

sc.stop()

