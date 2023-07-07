Novogene format for Sequencing Data

LINUX
Install Space Ranger 2.1.0 (May16 2023)

curl -o spaceranger-2.1.0.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-2.1.0.tar.gz?Expires=1686162097&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvc3BhdGlhbC1leHAvc3BhY2VyYW5nZXItMi4xLjAudGFyLmd6IiwiQ29uZGl0aW9uIjp7IkRhdGVMZXNzVGhhbiI6eyJBV1M6RXBvY2hUaW1lIjoxNjg2MTYyMDk3fX19XX0_&Signature=CoVVlUyToDTkHqatKKjhp6pERGUwhJhwqD1av1XB51CvuDTgOGeJqvv8Gl9RNQs-0hD~yEzUtmQLAJdYLBE8SnL-4rBybQz4ljp5W3cW5av1ggvDMVBeCmu3awQNKGZykcFDW4KzwKmVlSHqw9QoZcWE6gD7srOaeBTOOonYf~qJILRQLYGUSeNb9XmTq~8zulx-V1qx6Rb6BdkAsR~Qh6KwUbEybM-cdouxaQblTkAeNudZcRVTrsj87EVkGSi9x9paf5s7AmOlLXc3FsAiPHBgnkcq-RLKv3eZnivasBZu6Rsx1hxnDvDd4HxUf1CMy7b61T-2MJROj6STAvprvQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

Reference for human GRCh38 dataset required for Space Ranger

curl -O "https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-GRCh38-2020-A.tar.gz"

Download data from Novogene

wget -i ./X204SC23033668-Z01-F002.csv "https://objectstorage.uk-london-1.oraclecloud.com/p/8vuNRKi9O3cf91kvrRXd4U2KVzhy-8ZfM3juPo86hHWob-RKYkSjI3pXCwZXOVcx/n/cnyr09qj8zbo/b/england-data/o/out/CP2023030600059/X204SC23033668-Z01-F002/MD5.txt"

wget -i ./X204SC23033668-Z01-F002.csv "https://objectstorage.uk-london-1.oraclecloud.com/p/8vuNRKi9O3cf91kvrRXd4U2KVzhy-8ZfM3juPo86hHWob-RKYkSjI3pXCwZXOVcx/n/cnyr09qj8zbo/b/england-data/o/out/CP2023030600059/X204SC23033668-Z01-F002/X204SC23033668-Z01-F002.tar"

wget -i ./X204SC23033668-Z01-F002.csv "https://objectstorage.uk-london-1.oraclecloud.com/p/8vuNRKi9O3cf91kvrRXd4U2KVzhy-8ZfM3juPo86hHWob-RKYkSjI3pXCwZXOVcx/n/cnyr09qj8zbo/b/england-data/o/out/CP2023030600059/X204SC23033668-Z01-F002/checkSize.xls"

Install Space Ranger

https://www.10xgenomics.com/support/software/space-ranger/downloads/space-ranger-installation
Step 1 – Download and unpack the Space Ranger file in any location. This example uses /opt

$ cd /opt
[ download file from downloads page ]
$ tar -xzvf spaceranger-2.1.0.tar.gz

tar command
tar [options] [archive-file] [file or directory to be archived]
>>>> tar xzvf spaceranger-2.0.1.tar.gz 

Step2 - Download and unpack the reference of human Chr38

[ download file from downloads page ]
$ tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz

Step 3 – Prepend the Space Ranger directory to your $PATH. This allows you to invoke the spaceranger command.
If you unpacked both Space Ranger and the reference data into /opt, then run the following command as below. However, take a note that pwd should be inside the folder: spaceranger

$ vi ~/.bashrc

Copy paste this command to the file
$ export PATH=/home/leh/test/opt/spaceranger-2.1.0:$PATH

Working in Narvi
@ export PATH=/home/leh/test/opt/spaceranger-2.1.0:$PATH
After that, saving

$ source ~/.bashrc

Site check script
spaceranger upload hien.le@tuni.fi sitecheck.txt

Verify installation
To ensure that the spaceranger pipeline is installed correctly, take this test:
  
  spaceranger testrun --id=tiny


https://www.10xgenomics.com/support/software/space-ranger/analysis/inputs/input-overview#overview
Tracking down your inputs

scp -i ~/.ssh/narvi_key ./X204SC23033668-Z01-F002.tar leh@narvi.tut.fi:home/leh/data
The frist key to acess Narvi
ssh -i ~/.ssh/narvi_key leh@narvi.tut.fi

scp -i ~/.ssh/narvi_key ./X204SC23033668-Z01-F002.tar leh@narvi.tut.fi:/home/leh/data
scp -i ~/.ssh/narvi_key ./X204SC23033668-Z01-F002.tar leh@narvi.tut.fi:/home/leh/test

Spaceranger count
spaceranger count [FLAG] [OPTIONS] --id=sample1 --transcriptome=PATH --fastqs=PATH...
--image=IMG | --darkimage=IMG | --colorizedimage=IMG | --cytaimage=IMG
spaceranger count -h | --help | --version

unzip tar file in linux
tar -xf X204SC23033668-Z01-F002.tar

Upload images of capture areas to Narvi

scp -i ~/.ssh/narvi_key /home/hana/images_capture_area/capture_C1.png leh@narvi.tut.fi:/home/leh/test/images_capture_area
scp -i ~/.ssh/narvi_key /home/hana/images_capture_area/capture_C1.png leh@narvi.tut.fi:/home/leh/test/images_capture_area
upload new image: crop-C1

scp -i ~/.ssh/narvi_key /home/hana/images_capture_area/crop-C1.tiff leh@narvi.tut.fi:/home/leh/test/images_capture_area
scp -i ~/.ssh/narvi_key /home/hana/images_capture_area/crop-C1.tiff leh@narvi.tut.fi:/home/leh/test/images_capture_area

\\wsl.localhost\Ubuntu-22.04\home\hana\images_capture_area
spaceranger mkfastq
spaceranger mkfastq --run=PATH 

Automatic Alignment with spaceranger count

spaceranger count --ssh -i ~/.ssh/narvi_key leh@narvi.tut.fi
--transcriptome=/home/jdoe/refdata/GRCh38-2020-A \ #Path to Reference
--probe-set=/home/jdoe/spaceranger-2.0.0/probe_set/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv \ #Path to probe set
--fastqs=/home/jdoe/runs/HAWT7ADXX/outs/fastq_path \ #Path to FASTQs
--sample=mysample \ #Sample name from FASTQ filename
--image=/home/jdoe/runs/images/sample345.tiff \ #Path to brightfield image
--slide=V19J01-123 \ #Slide ID
--area=A1 \ #Capture area
--localcores=8 \ #Allowed cores in localmode
--localmem=64  #Allowed memory (GB) in localmode

cd /home/jdoe/runs
spaceranger count --id=/home/leh/test/capture_area_3
--transcriptome=/home/test/opt/refdata-gex-GRCh38-2020-A
--probe-set=/home/leh/test/opt/spaceranger-2.1.0/probe_set/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv                   --fastqs=/home/jdoe/runs/HAWT7ADXX/outs/fastq_path \ #Path to FASTQs
--sample=mysample \ #Sample name from FASTQ filename
--image=/home/jdoe/runs/images/sample345.tiff \ #Path to brightfield image
--slide=V11D13-338
--area=C1
--localcores=8
--localmem=64 


cd /home/leh
spaceranger count --id=/home/leh/test/capture_area_3 --transcriptome=/home/test/opt/refdata-gex-GRCh38-2020-A --fastqs=/home/leh/test/X204SC23033668-Z01-F002/01.RawData/Capture_area_3_3 --sample=Capture_area_3_3-AK4009-AK40100_HNNLKDSXS_S4_L003 --image=/home/jdoe/runs/images/sample345.tiff \ #Path to brightfield image --slide=V11D13-338 --area=C1 --localcores=8 --localmem=64 

spaceranger count --id capture_area_3 \
--transcriptome /home/leh//test/opt/refdata-gex-GRCh38-2020-A --image /home/leh/test/images_capture_area/capture_C1.png --fastqs /home/leh/test/X204SC23033668-Z01-F002/01.RawData/Capture_area_3_3/Capture_area_3_3-AK40099-AK40100_HNNLKDSX5_S4_L003_R1_001.fastq --sample Capture_area_3_3-AK4009-AK40100_HNNLKDSXS_S4_L003  --slide V11D13-338 --area C1 --localcores 8 --localmem 64 



spaceranger count pipeline will count gene expression and Feature Barcoding reads from a single capture area on the Visium slide.

Usage:
  
  scp -i ~/.ssh/narvi_key /home/hana/images_capture_area/capture_C1.png leh@narvi.tut.fi:/home/leh/test/images_capture_area

unzip file
gunzip Capture_area_3_3-AK40099-AK40100_HNNLKDSX5_S4_L003_R1_001.fastq.gz
gunzip Capture_area_3_3-AK40099-AK40100_HNNLKDSX5_S4_L003_R2_001.fastq.gz

22.6.2023: code can not run due to the dimension of images
ANH
spaceranger count --id capture_area_3 --transcriptome /home/leh/test/opt/refdata-gex-GRCh38-2020-A --image /home/leh/test/images_capture_area/capture_C1.png --fastqs /home/leh/test/X204SC23033668-Z01-F002/01.RawData/Capture_area_3_3 --sample Capture_area_3_3-AK40099-AK40100_HNNLKDSX5 --slide V11D13-338 --area C1 --localcores 8 --localmem 64



Result:[leh@narvi Capture_area_3_3]$ spaceranger count --id capture_area_3 --transcriptome /home/leh/test/opt/refdata-gex-GRCh38-2020-A --image /home/leh/test/images_capture_area/capture_C1.png --fastqs /home/leh/test/X204SC23033668-Z01-F002/01.RawData/Capture_area_3_3 --sample Capture_area_3_3-AK40099-AK40100_HNNLKDSX5 --slide V11D13-338 --area C1 --localcores 8 --localmem 64


Martian Runtime - v4.0.10
Serving UI at http://narvi.cc.tut.fi:43880?auth=oe5tWq3vKQ1Fd2H2Cu5KXVP6zDnEaIng9uwj56z1MGo

Running preflight checks (please wait)...
Checking sample info...
Checking FASTQ folder...
Checking reference...
Checking reference_path (/home/leh/test/opt/refdata-gex-GRCh38-2020-A) on narvi.cc.tut.fi...
Checking optional arguments...
mro: v4.0.10
mrp: v4.0.10
Anaconda: Python 3.10.10
numpy: 1.23.5
scipy: 1.10.1

[error] Image must have at least one dimension >= 2000 for standard slides

2023-06-23 08:32:18 Shutting down.
2023-06-23 08:32:18 Caught signal terminated
2023-06-23 08:32:18 [monitor] Caught signal terminated
Saving pipestance info to "capture_area_3/capture_area_3.mri.tgz"
For assistance, upload this file to 10x Genomics by running:
  
  spaceranger upload <your_email> "capture_area_3/capture_area_3.mri.tgz"

26.6.2023
spaceranger count --id capture_area_3 --transcriptome /home/leh/test/opt/refdata-gex-GRCh38-2020-A --image /home/leh/test/images_capture_area/crop-C1.tiff --fastqs /home/leh/test/X204SC23033668-Z01-F002/01.RawData/Capture_area_3_3 --sample Capture_area_3_3-AK40099-AK40100_HNNLKDSX5 --slide V11D13-338 --area C1 --localcores 8 --localmem 64

Upload image to narvi
scp -i ~/.ssh/narvi_key /home/hana/images_capture_area/capture_B1.tiff leh@narvi.tut.fi:/home/leh/test/images_capture_area

spaceranger count --id capture_area_2 --transcriptome /home/leh/test/opt/refdata-gex-GRCh38-2020-A --image /home/leh/test/images_capture_area/capture_B1.tiff --fastqs /home/leh/test/X204SC23033668-Z01-F002/01.RawData/Capture_area_2_2 --sample Capture_area_2_2-AK40097-AK40098_H2N3VDSX7 --slide V11D13-338 --area B1 --localcores 8 --localmem 64

gunzip Capture_area_2-2-AK40097-AK40098_H2N3VDSX7_S3_L001_I1_001.fastq.gz
gunzip Capture_area_2_2-AK40097-AK40098_H2N3VDSX7_S3_L001_I2_001.fastq.gz
gunzip Capture_area_2_2-AK40097-AK40098_H2N3VDSX7_S3_L001_R1_001.fastq.gz
gunzip Capture_area_2_2-AK40097-AK40098_H2N3VDSX7_S3_L001_R2_001.fastq.gz

27.6.2023 download data from serve to local computer
scp -i ~/.ssh/narvi_key leh@narvi.tut.fi:/home/leh/test/capture_area_2.tar.gz /home/hana/data
scp -i ~/.ssh/narvi_key leh@narvi.tut.fi:/home/leh/test/capture_area_3.tar.gz /home/hana/data
scp -i ~/.ssh/narvi_key leh@narvi.tut.fi:/home/leh/test/capture_area_1.tar.gz /home/hana/data

unzip tar file in linux
tar -xf capture_area_3.tar
tar -xf capture_area_2.tar
tar -xf capture_area_1.tar

Upload image to narvi
scp -i ~/.ssh/narvi_key /home/hana/images_capture_area/capture_A1.tiff leh@narvi.tut.fi:/home/leh/test/images_capture_area
scp -i ~/.ssh/narvi_key /home/hana/images_capture_area/capture_D1.tiff leh@narvi.tut.fi:/home/leh/test/images_capture_area

gunzip Capture_area_1_1-AK40095-AK40096_H2N3VDSX7_S2_L004_I1_001.fastq.gz
gunzip Capture_area_1_1-AK40095-AK40096_H2N3VDSX7_S2_L004_I2_001.fastq.gz
gunzip Capture_area_1_1-AK40095-AK40096_H2N3VDSX7_S2_L004_R1_001.fastq.gz
gunzip Capture_area_1_1-AK40095-AK40096_H2N3VDSX7_S2_L004_R2_001.fastq.gz
spaceranger count --id capture_area_1 --transcriptome /home/leh/test/opt/refdata-gex-GRCh38-2020-A --image /home/leh/test/images_capture_area/capture_A1.tiff --fastqs /home/leh/test/X204SC23033668-Z01-F002/01.RawData/Capture_area_1_1 --sample Capture_area_1_1-AK40095-AK40096_H2N3VDSX7 --slide V11D13-338 --area A1 --localcores 8 --localmem 64

gunzip Capture_area_4_4-AK40101-AK40102_HNNLKDSX5_S5_L004_I1_001.fastq.gz
gunzip Capture_area_4_4-AK40101-AK40102_HNNLKDSX5_S5_L004_I2_001.fastq.gz
gunzip Capture_area_4_4-AK40101-AK40102_HNNLKDSX5_S5_L004_R1_001.fastq.gz
gunzip Capture_area_4_4-AK40101-AK40102_HNNLKDSX5_S5_L004_R2_001.fastq.gz
spaceranger count --id capture_area_4 --transcriptome /home/leh/test/opt/refdata-gex-GRCh38-2020-A --image /home/leh/test/images_capture_area/capture_D1.tiff --fastqs /home/leh/test/X204SC23033668-Z01-F002/01.RawData/Capture_area_4_4 --sample Capture_area_4_4-AK40101-AK40102_HNNLKDSX5 --slide V11D13-338 --area D1 --localcores 8 --localmem 64

28.06.2023 
Download the data from narvi to local computer
scp -i ~/.ssh/narvi_key leh@narvi.tut.fi:/home/leh/test/capture_area_1 /home/hana/data (not sucess)

scp -i ~/.ssh/narvi_key leh@narvi.tut.fi:/home/leh/test/capture_area_1/outs/cloupe.cloupe /home/hana/data/capture_1
scp -i ~/.ssh/narvi_key leh@narvi.tut.fi:/home/leh/test/capture_area_1/outs/web_summary.html /home/hana/data/capture_1

scp -i ~/.ssh/narvi_key leh@narvi.tut.fi:/home/leh/test/capture_area_4/outs/cloupe.cloupe /home/hana/data/capture_4
scp -i ~/.ssh/narvi_key leh@narvi.tut.fi:/home/leh/test/capture_area_4/outs/web_summary.html /home/hana/data/capture_4

29.6.2023
unzip file in local computer
gunzip barcodes.tsv.gz
gunzip features.tsv.gz

2.7.2023 Working with R studio
data_cap_3 = Read10X_h5(filename = "filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
Load10X_Spatial(
  data.dir = "C:/Users/leh/OneDrive - TUNI.fi/Documents/Data/cap_3/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL,
)

3.7.2023 data processing following brain data in Seurat
plot1 <- VlnPlot(data_cap3, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(data_cap3, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

sctransform normalizes the data, detects high-variance features and stores the data in the SCT assay

data_cap3_sct = SCTransform(data_cap3, assay = "Spatial", verbose = FALSE)

Gene expression visualization
SpatialFeaturePlot(data_cap3_sct, features = c("GBP1"))

4.7.2023 Working in R (not R studio)

png(file="plot1.png")