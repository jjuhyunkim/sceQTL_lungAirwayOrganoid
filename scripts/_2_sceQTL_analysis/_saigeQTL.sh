celltype=$1
phenoCol=$2
export OMP_NUM_THREADS=1    # Limits OpenMP threads to 1
export MKL_NUM_THREADS=1    # Limits Intel MKL threads to 1 (used by NumPy, SciPy)
export OPENBLAS_NUM_THREADS=1  # Limits OpenBLAS threads to 1 (used by NumPy)
export NUMEXPR_NUM_THREADS=1   # Limits NumExpr threads to 1 (used by NumPy)
export VECLIB_MAXIMUM_THREADS=1  # Limits threads in macOS Accelerate framework (if applicable)


mainOutFolder=$3
mkdir -p $mainOutFolder/${celltype}/$phenoCol

pack=/data/User/juhyunk/project/02.sceQTL/saigeqtl/qtl/extdata/


wd=/gmi-l1/_90.User_Data/juhyunk/project/02.sceQTL/saigeqtl/run_lung_org
pheno=input/10.NRFonly.n11.merged.normMat.${celltype}.filter_genes.tsv


step1Plink=input/noLCR.maf005.lungNRF.n68.NARD2.hg38.mac20.ldprune.prune.in.radom1000
sampleid=sampleid
sampleCovar=gPC_1,gPC_2,gPC_3,gPC_4,gPC_5,sex
cov=ePC_3,ePC_4,ePC_5,log_nCount_RNA,$sampleCovar
threads=1
plinkFile=input/noLCR.maf005.lungNRF.n68.NARD2.hg38

step1out=$mainOutFolder/${celltype}/$phenoCol/lung_org_${celltype}_${phenoCol}_step1
step2commonout=$mainOutFolder/${celltype}/$phenoCol/lung_org_${celltype}_${phenoCol}_step2_common.txt
step2rareout=$mainOutFolder/${celltype}/$phenoCol/lung_org_${celltype}_${phenoCol}_step2_rare.txt
step3commonout=$mainOutFolder/${celltype}/$phenoCol/lung_org_${celltype}_${phenoCol}_step3_common.txt

# step1 
if [ ! -f ${step1out}.rda ]; then
echo "run step1 : ${step1out}.rda"

Rscript $pack/step1_fitNULLGLMM_qtl.R \
--useSparseGRMtoFitNULL=False  \
--useGRMtoFitNULL=False \
--phenoFile=${pheno} \
--phenoCol=$phenoCol \
--covarColList=$cov \
--sampleCovarColList=$sampleCovar \
--sampleIDColinphenoFile=$sampleid \
--traitType=count \
--outputPrefix=${step1out} \
--skipVarianceRatioEstimation=FALSE \
--isRemoveZerosinPheno=FALSE \
--isCovariateOffset=FALSE \
--isCovariateTransform=FALSE \
--skipModelFitting=FALSE \
--tol=0.00001 \
--plinkFile=$step1Plink \
--maxiterPCG=10000 \
--nThreads=$threads \
--IsOverwriteVarianceRatioFile=TRUE
else
echo "skip step1"
fi

# step2 common cis
regionFile=input/regions/$phenoCol.region.txt

if [ ! -f ${regionFile} ]; then
echo "RegionFile generation"
# fgrep -w $phenoCol input/10.NRFonly.n11.merged.normMat.filter_genes.region_padd_1Mbp.tsv  | cut -f 1-3 > $regionFile 
fgrep -w $phenoCol /data/User/juhyunk/project/02.sceQTL/saigeqtl/run_lung_org/input/hLRO_merged.3DNetMod10k.domain.filtered_blacklist.byGene.bed | cut -f 1-3 | sort | uniq  | sed -e 's/^chr//' > $regionFile
fi

if [ ! -f ${step2commonout} ]; then
echo "step2 common"
Rscript $pack/step2_tests_qtl.R \
--bedFile=$plinkFile.bed \
--bimFile=$plinkFile.bim \
--famFile=$plinkFile.fam \
--SAIGEOutputFile=${step2commonout} \
--minMAF=0 \
--minMAC=1 \
--LOCO=FALSE \
--GMMATmodelFile=${step1out}.rda \
--SPAcutoff=2 \
--varianceRatioFile=${step1out}.varianceRatio.txt \
--markers_per_chunk=10000 \
--rangestoIncludeFile=$regionFile
fi

# step2 rare set based 


groupFile=${regionFile}.grp

if [ ! -f ${groupFile} ]; then
echo "groupFile generation"
Rscript $pack/makeGroupFile.R \
--bedFile=$plinkFile.bed \
--bimFile=$plinkFile.bim \
--famFile=$plinkFile.fam \
--regionFile=${regionFile} \
--outputPrefix=$groupFile 
fi

if [ ! -f ${step2rareout} ]; then
echo "$step2rareout"
Rscript $pack/step2_tests_qtl.R \
--bedFile=$plinkFile.bed \
--bimFile=$plinkFile.bim \
--famFile=$plinkFile.fam \
--SAIGEOutputFile=${step2rareout} \
--maxMAF_in_groupTest=0.1 \
--minMAF_in_groupTest_Exclude=0 \
--groupFile=${groupFile} \
--annotation_in_groupTest=null \
--MACCutoff_to_CollapseUltraRare=10 \
--is_single_in_groupTest=TRUE \
--is_equal_weight_in_groupTest=TRUE \
--LOCO=FALSE \
--GMMATmodelFile=${step1out}.rda \
--SPAcutoff=2 \
--varianceRatioFile=${step1out}.varianceRatio.txt \
--markers_per_chunk=10000
fi

# step3 gene based
if [ ! -f ${step3commonout} ]; then
echo "${step3commonout}"
Rscript $pack/step3_gene_pvalue_qtl.R \
        --assocFile=${step2commonout} \
        --geneName=$phenoCol \
        --genePval_outputFile=$step3commonout
fi
