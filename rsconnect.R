# Deploy to shinyapps.io. Credentials are read from environment variables:
#   SHINYAPPS_ACCOUNT, SHINYAPPS_TOKEN, SHINYAPPS_SECRET
# e.g. put them in ~/.Renviron, never in this file.
library(rsconnect)

account <- Sys.getenv("SHINYAPPS_ACCOUNT", "crafor18")
token <- Sys.getenv("SHINYAPPS_TOKEN")
secret <- Sys.getenv("SHINYAPPS_SECRET")
if (!nzchar(token) || !nzchar(secret)) {
  stop("Set SHINYAPPS_TOKEN and SHINYAPPS_SECRET before deploying.")
}

rsconnect::setAccountInfo(name = account, token = token, secret = secret)
rsconnect::deployApp(
  appDir = ".",
  appName = "Z_Y_ketang_data_analysis_app",
  appFiles = c("app.R", "R", "www", "docs", "data/examples", "fonts")
)
