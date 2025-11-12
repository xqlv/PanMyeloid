seurat_path=~/seurat_obj.rds
output_path=~/cMO

seurat_path_basename=`basename $seurat_path`
prefix=${seurat_path_basename/.rds/}
output_path=${output_path}/${prefix}
mkdir -p $output_path

loom_out=${output_path}/expr_mat.loom
grn_out=${output_path}/expr_mat.adjacencies.tsv
ctx_out=${output_path}/regulons.csv
auc_out=${output_path}/auc_mtx.csv

Rscript ~/seurat2loom.R -i $seurat_path -o $loom_out

/opt/pyscenic/0.12.1/bin/pyscenic grn \
--num_workers 5 \
--seed 777 \
-o $grn_out \
$loom_out \
/DATA/public/cistarget/tf_lists/allTFs_hg38.txt

/opt/pyscenic/0.12.1/bin/pyscenic ctx \
$grn_out \
/DATA/public/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
/DATA/public/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather   \
--annotations_fname /DATA/public/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname $loom_out \
--expression_mtx_fname $loom_out \
--mode dask_multiprocessing \
--output $ctx_out \
--num_workers 5


/opt/pyscenic/0.12.1/bin/pyscenic aucell \
$loom_out \
$ctx_out \
-o $auc_out \
--num_workers 5

/opt/pyscenic/0.12.1/bin/python ~/regulons2df.py \
$ctx_out \
$output_path/tfs_targets.csv 
