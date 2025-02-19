# 安装 rsconnect 包（如未安装）
# install.packages("rsconnect")
library(rsconnect)
setwd("D:\\曾昀毕业论文\\8-大论文\\GWAS论文修改\\从头开始的分析\\tools\\merge_shinny")
# 设置账号信息（首次部署时需要配置）
rsconnect::setAccountInfo(name='crafor18',
                          token='776199A3E6A0AF06EDD3B33982ACA2EB',
                          secret='DccvEttJD0r0X4WmuWJJNfZ+beea6W+tpu8FYRRd')

# 部署应用
rsconnect::deployApp(appDir = ".", appName = "Z_Y_ketang_data_analysis_app")
