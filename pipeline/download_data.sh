#!/bin/bash

##### ---------- inhouse data (from huaweicloud) ---------- #####
./obsutil share-cp https://e-share.obs-website.cn-north-1.myhuaweicloud.com?token=xxxxxxxxxx公司返回数据的邮件中很长的授权码（华为云文件夹链接） xxxxxxx/xxxxxxxx/服务器上储存数据的目标文件夹/ -r -vmd5


##### ---------- downloading public data ---------- #####
cd $workdir
download_raw_seq.py -p $PRJNAXXXXXX -f 'SRRxxxxxx1|SRRxxxxxx2|SRRxxxxxx3' -o 0_raw_data
cd ${workdir}/0_raw_data
bash ${PRJNAXXXXXX}.sh