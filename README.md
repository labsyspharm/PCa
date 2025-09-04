# Data Access Information
------------------ 
TABLE OF CONTENTS
------------------
​
* GENERAL INFORMATION
  * ASSOCIATED PUBLICATION
  * RECOMMENDED CITATION
  * USEFUL LINKS
* ACCESS THE DATASET
  * FILE ORGANIZATION
  * REPOSITORY LINKS
  * FILE LIST
* ADDITIONAL NOTES/COMMENTS
​
--------------------
GENERAL INFORMATION
--------------------
​
1. **Spatial organization of the lymphocyte infiltration in prstate tumors** <Publication or Dataset Title>   
​
2. **Authors:** <List; No need to include affiliations.>  
​
3. **Please cite this data as the following:**  <Following standard APA citation format - if this data is associated with a publication, this should be the paper citation>    
Author Last, Author F. (Year). Title of data set (Version number) [Description of form]. Location: Name of producer.    
  
4. **Relevant links:** <remove links that are not relevant>  
> * Publication DOI: [doi.org/MY-PAPER-DOI](https://doi.org/MY-PAPER-DOI-URL) 
> * Associated GitHub Repository: [MY-REPO](https://github.com/MY-REPO-URL)  
> * To view an archived record of this repository: [My-ZENODO-DOI](https://zenodo.org/doi/MY-ZENDODO-DOI-URL) 
> * To view the image data online, visit: [My-ATLAS-PAGE](https://tissue-atlas.org/MY-ATLAS-PAGE-URL)
> * Associated documentation website: [MY-PROJECT-WEBSITE](MY-PROJECT-WEBSITE-URL)  
> * DOI of other publications that use the data: <If this data is being reused from a past publication, include DOI and APA citation>
> * <Other important links if applicable>
​
5. **Licenses/restrictions placed on the data:** CC-BY [creativecommons.org/licenses/by/4.0/](https://creativecommons.org/licenses/by/4.0/)
​
--------------------
ACCESS THE DATASET 
--------------------
<Data should be uploaded to the appropriate public repository where applicable - if you are not sure which repository to use reach out to your Data Manager>
  
​
## File organization:   
**Each file follows the following naming convention:**    
​
Each folder corresponds to a patient sample (N). <Edit as needed if this folder structure does not fit the needs of your paper> 
 
|File Type     | Description                                                                        | Location|
|--------      | ----------------------------------------------------------------------------------|---------|
|N.ome.tif     | Stitched multiplex CyCIF image pyramid in ome.tif format                           | AWS     |
|N_HE.vsi      | Hematoxylin and Eosin stained image of adjacent FFPE tissue section in .vsi format | AWS     |
|\_N\_HE\_/    | folder: raw image data accompanying .vsi file                                      | AWS     |
|markers.csv   | list of all markers in ome.tif image                                               | Synapse |
|N.csv         | single-cell feature table, including intensity data for all channels               | Synapse |
|N_ROI.csv     | X and Y coordinates for histologically annotated regions in CyCIF and H&E images   | Synapse |
|raw/          | folder of raw IF image data in .rcpnl format                                       | AWS     |
|segmentation/ |  folder of segmentation maps for tissue image in .tif format                       | AWS     |
| N.fastq      | Sequencing data                                                                    | GEO     |
​
​
## AWS Data Access  
​
X, Y, Z data is available for download through AWS. 
​
**You will need the following bucket name:**  
```
AWS BUCKET NAME  
```
​
*For general instructions on how to download data from AWS, see: [https://zenodo.org/records/10223574](https://zenodo.org/records/10223574)*     
  
If you experience issues accessing the above AWS S3 buckets, email tissue-atlas(at)hms.harvard.edu with the subject line "bucketname: Data Access".  
​
## Synapse Library
​
X, Y, Z, data is available on Synapse at: [synapse.org/MY-PUBLIC-SYNAPSE-LIBRARY](https://www.synapse.org/MY-PUBLIC-SYNAPSE-LIBRARY)  
​
## Other Repository 
<if applicable>  
​
## FILE LIST  
<List all files (or folders, as appropriate for dataset organization) contained in each repository, with a brief description. If you are depositing certain file types into public, standardized repositories that already include a file index & metadata, you can link to that repository instead of listing all individual files. For all other data, (on AWS, etc) list all files. >  
​
### N.ome.tif
​
|Patient or Biospecimen ID | File Name       | Location| File size |
|------- | ----------------|---------|-----------|
|ID | ID.ome.tif | AWS     | N.N GB   |
​
​
### N_HE.vsi
​
|Patient or Biospecimen ID | File Name      | Location| File size|
|--------| ---------------|---------|----------|
|ID | ID.ome.tif_HE.vsi | AWS     | N.N KB |
​
​
### \_N\_HE\_/
​
|Patient or Biospecimen ID | File Name   | Location| File size|
|------- | ------------|---------|----------|
|ID | frame_t.ets | AWS     | N.N MB |
​
### markers.csv
​
|Patient or Biospecimen ID | File Name   | Synapse ID  | File size|
|------- | ----------- |------------ |----------|
|ID | markers.csv | syn12345678 | N.N bytes|
​
### N.csv
​
|Patient or Biospecimen ID | File Name   | Synapse ID | File size |
|------- | ------------|------------|-----------|
|ID | ID.csv |  |  |
​
### N_ROI.csv
​
|Patient or Biospecimen ID | File Name       | Synapse ID  | File size|
|------- | ----------------|-------------|----------|
|ID | ID_ROI.csv |  |    |
​
### raw/
​
|Patient or Biospecimen ID | Number of Files | Folder size| Location|
|------- |-----------------|------------|---------|
|ID | 13              |     |      |
​
​
### segmentation/
​
|Patient or Biospecimen ID | Number of Files | Folder size| Location|
|------- |-----------------|------------|---------|
|ID |               |     |      |
​
​
 
--------------------------
ADDITIONAL NOTES/COMMENTS
--------------------------
​
Please let **(Name, contact information)** know if any errors are found in this data.  
