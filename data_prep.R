# SETUP ----
if(!"easypackages" %in% installed.packages()[, "Package"]) {
	install.packages("easypackages")
}

library(easypackages)
packages("tidyverse", "tidyquant")

# This is a helper file, to be used simply to prepare the data for the shiny app

# Data Download ---
from <- Sys.Date() - years(40)

stocks <- tq_index("SP500", use_fallback = T)["symbol"] %>% 
	tq_get(get = "stock.prices", from = from) %>% 
	group_by(symbol) %>% 
	distinct(date, .keep_all = T) %>% 
	ungroup()

snp <- tq_get("^GSPC", get = "stock.prices", from = from) %>% 
	distinct(date, .keep_all = T)

# 10 year yield is conventionally used as the risk-free rate
rf <- tq_get("DGS10", get = "economic.data", from = from + weeks(1)) %>% 
	distinct(date, .keep_all = T)

# Data Wrangling ----
# Set data frequencies
d_rets <- stocks %>% 
	dplyr::select(symbol, date, adjusted) %>% 
	group_by(symbol) %>% 
	transmute(
		date = date,
		ret = adjusted / lag(adjusted) - 1,
		freq = "daily"
	) %>% 
	drop_na() %>% 
	filter(n() > 1) %>% 
	ungroup()

d_snp_ret <- snp %>% 
	select(date, adjusted) %>% 
	transmute(
		symbol = "GSPC",
		date = date,
		ret = adjusted / lag(adjusted) - 1,
		freq = "daily"
	) %>% 
	drop_na()

d_rf_ret <- rf %>% 
	# Yield is annualized, so it is divided it accordingly to make it comparable
	# to the other return series; it is also divided by 100, as the raw data is
	# provided in percentage units
	transmute(
		symbol = "DGS10",
		date = date,
		ret = price / 100 / 365,
		freq = "daily"
	) %>% 
	drop_na()

w_rets <- stocks %>% 
	dplyr::select(symbol, date, adjusted) %>% 
	group_by(symbol) %>% 
	tq_transmute(
		select = adjusted,
		mutate_fun = to.weekly,
		indexAt = "lastof"
	) %>%
	transmute(
		date = date,
		ret = adjusted / lag(adjusted) - 1,
		freq = "weekly"
	) %>% 
	drop_na() %>% 
	filter(n() > 1) %>% 
	ungroup()

w_snp_ret <- snp %>% 
	select(date, adjusted) %>% 
	tq_transmute(
		select = adjusted,
		mutate_fun = to.weekly,
		indexAt = "lastof"
	) %>% 
	transmute(
		symbol = "GSPC",
		date = date,
		ret = adjusted / lag(adjusted) - 1,
		freq = "weekly"
	) %>% 
	drop_na()

w_rf_ret <- rf %>% 
	tq_transmute(
		select = price,
		mutate_fun = to.weekly,
		indexAt = "lastof"
	) %>%
	transmute(
		symbol = "DGS10",
		date = date,
		ret = price / 100 / 52,
		freq = "weekly"
	) %>% 
	drop_na()

m_rets <- stocks %>% 
	dplyr::select(symbol, date, adjusted) %>% 
	group_by(symbol) %>% 
	tq_transmute(
		select = adjusted,
		mutate_fun = to.monthly,
		indexAt = "lastof"
	) %>%
	transmute(
		date = date,
		ret = adjusted / lag(adjusted) - 1,
		freq = "monthly"
	) %>% 
	drop_na() %>% 
	filter(n() > 1) %>% 
	ungroup()

m_snp_ret <- snp %>% 
	select(date, adjusted) %>% 
	tq_transmute(
		select = adjusted,
		mutate_fun = to.monthly,
		indexAt = "lastof"
	) %>% 
	transmute(
		symbol = "GSPC",
		date = date,
		ret = adjusted / lag(adjusted) - 1,
		freq = "monthly"
	) %>% 
	drop_na()

m_rf_ret <- rf %>% 
	tq_transmute(
		select = price,
		mutate_fun = to.monthly,
		indexAt = "lastof"
	) %>%
	transmute(
		symbol = "DGS10",
		date = date,
		ret = price / 100 / 12,
		freq = "monthly"
	) %>% 
	drop_na()

q_rets <- stocks %>% 
	dplyr::select(symbol, date, adjusted) %>% 
	group_by(symbol) %>% 
	tq_transmute(
		select = adjusted,
		mutate_fun = to.quarterly,
		indexAt = "lastof"
	) %>%
	transmute(
		date = date,
		ret = adjusted / lag(adjusted) - 1,
		freq = "quarterly"
	) %>% 
	drop_na() %>% 
	filter(n() > 1) %>% 
	ungroup()

q_snp_ret <- snp %>% 
	select(date, adjusted) %>% 
	tq_transmute(
		select = adjusted,
		mutate_fun = to.quarterly,
		indexAt = "lastof"
	) %>% 
	transmute(
		symbol = "GSPC",
		date = date,
		ret = adjusted / lag(adjusted) - 1,
		freq = "quarterly"
	) %>% 
	drop_na()

q_rf_ret <- rf %>% 
	tq_transmute(
		select = price,
		mutate_fun = to.quarterly,
		indexAt = "lastof"
	) %>%
	transmute(
		symbol = "DGS10",
		date = date,
		ret = price / 100 / 4,
		freq = "quarterly"
	) %>% 
	drop_na()

y_rets <- stocks %>% 
	dplyr::select(symbol, date, adjusted) %>% 
	group_by(symbol) %>% 
	tq_transmute(
		select = adjusted,
		mutate_fun = to.yearly,
		indexAt = "lastof"
	) %>%
	transmute(
		date = date,
		ret = adjusted / lag(adjusted) - 1,
		freq = "yearly"
	) %>% 
	drop_na() %>% 
	filter(n() > 1) %>% 
	ungroup()

y_snp_ret <- snp %>% 
	select(date, adjusted) %>% 
	tq_transmute(
		select = adjusted,
		mutate_fun = to.yearly,
		indexAt = "lastof"
	) %>% 
	transmute(
		symbol = "GSPC",
		date = date,
		ret = adjusted / lag(adjusted) - 1,
		freq = "yearly"
	) %>% 
	drop_na()

y_rf_ret <- rf %>% 
	tq_transmute(
		select = price,
		mutate_fun = to.yearly,
		indexAt = "lastof"
	) %>%
	transmute(
		symbol = "DGS10",
		date = date,
		ret = price / 100,
		freq = "yearly"
	) %>% 
	drop_na()

data <- bind_rows(
	d_rets, d_snp_ret, d_rf_ret,
	w_rets, w_snp_ret, w_rf_ret,
	m_rets, m_snp_ret, m_rf_ret,
	q_rets, q_snp_ret, q_rf_ret,
	y_rets, y_snp_ret, y_rf_ret
)

write_csv(data, "data.csv")

# Data for user download ----
stocks2 <- stocks %>% 
	select(symbol, date, adjusted) %>% 
	rename(obs = adjusted)

snp2 <- snp %>% 
	select(date, adjusted) %>% 
	add_column(symbol = "GSPC", .before = "date") %>% 
	rename(obs = adjusted)

rf2 <- rf %>% 
	transmute(
		date = date,
		obs = price / 100
	) %>% 
	add_column(symbol = "DGS10", .before = "date")

data2 <- bind_rows(stocks2, snp2, rf2)

write_csv(data2, "data2.csv")
