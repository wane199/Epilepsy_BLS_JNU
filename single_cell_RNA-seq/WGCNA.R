# https://mp.weixin.qq.com/s?__biz=MzAwMzIzOTk5OQ==&mid=2247499086&idx=1&sn=85f1a292ab82fc3c5980c82069057a79&chksm=9b3c841eac4b0d08bb8b821f67b2740a677b040d310e414c28e1fec7bd2a5cf4e0cef4ba0b0f&mpshare=1&scene=1&srcid=04225hR6eUh2LxTU03RZyBTz&sharer_sharetime=1682147335292&sharer_shareid=13c9050caaa8b93ff320bbf2c743f00b#rd
rm(list = ls())
getwd()
# 预处理
# 1.1 读取数据
if(!require(WGCNA))BiocManager::install("WGCNA")