# ##############################################################################
# SETUP
# ##############################################################################

library(tidyverse)
library(modelr)
library(broom)
library(tidyquant)
library(timetk)
library(tibbletime)
library(quadprog)
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(R6)
library(Matrix)

data <- read_csv("data.csv")
data2 <- read_csv("data2.csv")
min_date <- Sys.Date() - years(40)
max_date <- Sys.Date() - wday(Sys.Date() + 1) - days(7)

# Number of assets (for data wrangling purposes)
n_assets <- data %>% distinct(symbol) %>% nrow(.)

# Set the global names of optimization approaches
optim_names <- list("Single Index", "Constant Correlation",
								    "Minimum Variance with Target Return")
									  
optim_values <- list("SingleIndex", "ConsCorr", "MinVar")

eda_plot_heights <- 250

# ##############################################################################
# DESIGN OF R6 CLASSES
# ##############################################################################

SingleIndex <- 
	R6Class(
		"SingleIndex",
		private = list(
			# A function to model betas
			.fn_si = function(df) {
				lm(ret.i ~ ret.m, data = df, na.action = "na.exclude")
			},
			
			# A function to calculate portfolio variance
			.fn_p_var = function(rets_df, weights) {
				rets_df %>% 
					cov(., use = "pairwise.complete.obs") %>% 
					as.vector() %>% 
					`*`(., rep(weights, self$n_pos)) %>% 
					`*`(., rep(weights, each = self$n_pos)) %>% 
					sum(., na.rm = T)
			},
			
			# A function to calculate portfolio correlation
			.fn_p_rho = function(rets_df) {
				rets_df %>% 
					cor(., use = "pairwise.complete.obs") %>% 
					.[lower.tri(.)] %>% 
					mean(., na.rm = T)
			},
			
			.backtest_portfolios = NA,
			.backtest_output = NA,
			.aa_portfolio = NA,
			.aa_output = NA
		), # end of private list
		
		public = list(
			initialize = function(
				r_tibble = NULL, n_pos = NULL, roll_win = NULL, ...
			) {
				self$r_tibble <- r_tibble
				self$n_pos <- n_pos
				self$roll_win <- roll_win
			}, # end of initialize()
			
			r_tibble = NULL,
			n_pos = 10000,
			roll_win = 52,
			
			backtest = function() {
				
				# Set the rolling functions
				rolling_si <- rollify(.f = function(ret.i, ret.m) {
					lm(ret.i ~ ret.m)
				},
				window = self$roll_win,
				unlist = F
				)
				
				rolling_mean <- rollify(.f = ~mean(.x, na.rm = T),
																window = self$roll_win
				)
				
				rolling_sd <- rollify(.f = ~sd(.x, na.rm = T),
															window = self$roll_win
				)
				
				# Run the backtest
				rets <- self$r_tibble %>% 
					group_by(symbol) %>% 
					
					# For rolling functions, observations must be at least as big as the 
					# rolling window
					filter(n() >= self$roll_win) %>% 
					
					# Add residuals - necessary for 'Single Index' optimization
					nest() %>% 
					mutate(
						model = map(data, private$.fn_si),
						resids = map2(data, model, add_residuals)
					) %>% 
					select(symbol, resids) %>% 
					unnest(cols = resids) %>%
					
					# Calculate the statistics required for 'Single Index' optimization
					mutate(
						roll.si = rolling_si(ret.i, ret.m), # To get assets' beta
						exp.i = rolling_mean(ret.i), # Asset expected return
						std.e = rolling_sd(resid), # Asset nonsystematic risk
						std.m = rolling_sd(ret.m) # Market std dev
					) %>% 
					filter(!is.na(exp.i)) %>% # Remove unnecessary rows
					mutate(tidied = map(roll.si, tidy)) %>%
					unnest(tidied) %>%
					filter(term == "ret.m")%>% # keep only the beta term
					rename(beta = estimate) %>% 
					# For now, take out the p.value - will be necessary if the model is 
					# to be adjusted for high p-values
					select(
						date, symbol, ret.i, exp.i, beta, std.e, ret.m, std.m, ret.f
					) %>% 
					ungroup() %>% 
					arrange(date) %>% 
					group_by(date) %>% 
					
					# Compute ratio of excess return to beta, then rank them by 
					# desirability
					mutate(excess = (exp.i - ret.f) / beta) %>% 
					.[!is.infinite(.$excess), ] %>% 
					arrange(desc(excess), .by_group = T) %>% 
					
					# Compute the cutoff rate coefficient to then calculate c.star, the 
					# cutoff rate
					mutate(
						c.coef = (std.m^2 * cumsum((exp.i - ret.f) * beta / std.e^2)) /
							(1 + std.m^2 * cumsum(beta^2 / std.e^2))
					)
				
				# Extract the cutoff rate, the point at which the assets above it are 
				# expected to contribute to the excess returns of the optimal 
				# portfolio 
				c.star <- rets %>% 
					filter(excess >= c.coef, .preserve = T) %>% 
					top_n(1, c.coef) %>% 
					rename(c.star = c.coef) %>% 
					select(date, c.star)
				
				# Apply the cutoff rate
				rets <- rets %>% 
					left_join(c.star, by = "date") %>% 
					mutate(c.star = replace_na(c.star, 0),
								 
								 # Compute the investment proportion coefficient
								 z.coef = (beta * std.e^2) * (excess - c.star)
					) %>% 
					
					# Allowing for short sales, rank assets by profitability potential
					arrange(desc(abs(z.coef)), .by_group = T) %>% 
					
					# Apply the cardinality constraint
					top_n(self$n_pos, abs(z.coef)) %>% 
					
					# Apply Lintner's short sale definition
					mutate(
						weight = z.coef / sum(abs(z.coef)),
						p.beta = sum(beta * abs(weight))
					) %>% 
					arrange(symbol, .by_group = T) %>% 
					ungroup()
				
				# Compute portfolio variances and mean correlations
				# Get the symbols in each period's portfolio
				asset_names <- rets %>% 
					select(date, symbol) %>% 
					pivot_wider(names_from = symbol, values_from = symbol) %>% 
					ungroup() %>% 
					arrange(date) %>% 
					select(date, order(colnames(.)))
				
				# Format returns tibble
				rets_df <- self$r_tibble %>% 
					select(symbol, date, ret.i) %>% 
					pivot_wider(names_from = symbol, values_from = ret.i) %>% 
					select(date, order(colnames(.))) %>%
					ungroup() %>% 
					distinct(date, .keep_all = T) %>% 
					arrange(date)
				
				# Generate the rolling windows
				idx <- map(
					seq(1, nrow(rets_df) - self$roll_win, by = 1), 
					function(i) seq(i, i + self$roll_win - 1)
				)
				
				# Get the returns data frame that reflects each period's optimized 
				# portfolio
				roll_rets <- map(idx, function(i) rets_df[i, ])
				
				roll_rets <- map(
					seq_along(roll_rets),
					function(i)
						roll_rets[[i]][, colnames(roll_rets[[i]]) %in% asset_names[i, ]]
				)
				
				# Compute the variance-covariance matrices
				roll_cov <- map(
					seq_along(roll_rets),
					function(i)
						cov(roll_rets[[i]], use = "pairwise.complete.obs")
				)
				
				# Get corresponding weights
				asset_weights <- rets$weight %>% 
					matrix(., ncol = self$n_pos, byrow = T)
				
				# Calculate portfolio variance
				p_var <- map_dbl(
					seq_along(roll_cov),
					function(i) {
						sum(
							rep(asset_weights[i, ], self$n_pos) * 
								rep(asset_weights[i, ], each = self$n_pos) *
								as.vector(roll_cov[[i]])	
						)
					}
				)
				
				# Compute the correlation matrices
				roll_rho <- map(
					seq_along(roll_rets),
					function(i)
						cor(roll_rets[[i]], use = "pairwise.complete.obs")
				)
				
				# Calculate mean portfolio correlations
				p_rho <- map_dbl(
					seq_along(roll_rho),
					function(i) {
						roll_rho[[i]][lower.tri(roll_rho[[i]])] %>% mean(., na.rm = T)
					}
				)
				
				rets <- rets %>% 
					add_column(
						p.var = c(
							rep(p_var, each = self$n_pos),
							if(length(p_var) * self$n_pos < nrow(rets))
								rep(NA, self$n_pos)
						)
					) %>% 
					add_column(
						p.rho = c(
							rep(p_rho, each = self$n_pos),
							if(length(p_rho) * self$n_pos < nrow(rets))
								rep(NA, self$n_pos)
						)
					)
				
				# Compute portfolio performance and format appropriately
				portfolios <- rets %>% 
					group_by(symbol) %>% 
					
					# A rolling window backtest performance is computed with lagged weights
					mutate(
						ret.b = lag(weight, order_by = date) * ret.i,
						ret.b = replace_na(ret.b, 0), # Helpful for graphing purposes
						m.var = std.m^2
					) %>% 
					ungroup() %>% 
					group_by(date) %>% 
					mutate(ret.p = sum(ret.b)) %>% 
					
					# Format the portfolio tibble
					arrange(date) %>% 
					distinct(date, .keep_all = T) %>% 
					select(date, ret.p, p.beta, p.var, p.rho, ret.m, m.var) %>% 
					rename(SingleIndex = ret.p, index = ret.m) %>% 
					gather(index, SingleIndex, key = "set", value = "return") %>% 
					ungroup() %>% 
					group_by(set) %>% 
					
					# Compute compound returns
					mutate(
						return = replace(return, row_number() == 1, 0),
						level = cumprod(1 + return),
						p.beta = lag(p.beta),
						p.var = lag(p.var),
						p.rho = lag(p.rho)
					) %>% 
					ungroup()
				
				private$.backtest_portfolios  <- portfolios
				private$.backtest_output <- rets
			}, # end of backtest()
			
			asset_allocation = function() {
				
				# Get the last period risk-free rate
				rf_last <- self$r_tibble %>% 
					filter(date == max(date)) %>% 
					distinct(ret.f) %>% 
					as.numeric(.)
				
				# Run the asset allocation
				rets_tibble <- self$r_tibble %>% 
					select(symbol, date, ret.i) %>% 
					pivot_wider(names_from = symbol, values_from = ret.i) %>% 
					arrange(desc(date)) %>% 
					top_n(self$roll_win, wt = date) %>% 
					pivot_longer(2:ncol(.), names_to = "symbol", values_to = "ret.i") %>% 
					drop_na() %>% 
					group_by(symbol) %>% 
					arrange(date, .by_group = T) %>% 
					ungroup()
				
				to_join <- self$r_tibble %>% 
					select(date, ret.m, ret.f) %>% 
					distinct(date, .keep_all = T)
				
				rets <- rets_tibble %>% 
					group_by(symbol) %>% 
					left_join(to_join, by = "date") %>% 
					filter(n() > 1) %>% 
					
					# Add residuals - necessary for 'Single Index' optimization
					nest() %>% 
					mutate(
						model = map(data, private$.fn_si),
						resids = map2(data, model, add_residuals),
						tidied = map(model, tidy)
					) %>% 
					unnest(cols = tidied) %>%
					filter(term == "ret.m")%>% # keep only the beta term
					rename(beta = estimate) %>% 
					select(symbol, resids, beta) %>% 
					unnest(cols = resids) %>% 
					
					# Calculate the statistics required for 'Single Index' optimization
					mutate(
						exp.i = mean(ret.i, na.rm = T), # Asset expected return
						std.i = sd(ret.i, na.rm = T), # Asset standard deviation
						std.e = sd(resid, na.rm = T), # Asset nonsystematic risk
						std.m = sd(ret.m, na.rm = T), # Market std dev
						ret.f = rf_last
					) %>% 
					distinct(symbol, .keep_all = T) %>% 
					# For now, take out the p.value - will be necessary if the model is 
					# to be adjusted for high p-values
					select(
						symbol, ret.i, exp.i, beta, std.i, std.e, ret.m, std.m, ret.f
					) %>% 
					ungroup() %>% 
					
					# Compute ratio of excess return to beta, then rank them by 
					# desirability
					mutate(excess = (exp.i - ret.f) / beta) %>% 
					.[!is.infinite(.$excess), ] %>% 
					arrange(desc(excess), .by_group = T) %>% 
					
					# Compute the cutoff rate coefficient to then calculate c.star, the 
					# cutoff rate
					mutate(
						c.coef = (std.m^2 * cumsum((exp.i - ret.f) * beta / std.e^2)) /
							(1 + std.m^2 * cumsum(beta^2 / std.e^2))
					)
				
				# Extract the cutoff rate, the point at which the assets above it are 
				# expected to contribute to the excess returns of the optimal 
				# portfolio 
				c.star <- rets %>% 
					filter(excess >= c.coef, .preserve = T) %>% 
					top_n(1, c.coef) %>% 
					rename(c.star = c.coef) %>% 
					select(c.star)
				
				# Apply the cutoff rate
				rets <- rets %>% 
					add_column(c.star = c.star$c.star) %>% 
					
					# Compute the investment proportion coefficient
					mutate(z.coef = (beta * std.e^2) * (excess - c.star)) %>% 
					
					# Allowing for short sales, rank assets by profitability potential
					arrange(desc(abs(z.coef)), .by_group = T) %>% 
					
					# Apply the cardinality constraint
					top_n(self$n_pos, abs(z.coef)) %>% 
					
					# Apply Lintner's short sale definition
					mutate(
						weight = z.coef / sum(abs(z.coef)),
						exp.p = sum(weight * exp.i),
						p.beta = sum(beta * abs(weight))
					) %>% 
					arrange(symbol, .by_group = T) %>% 
					ungroup()
				
				# Compute portfolio variances and mean correlations
				# Get the symbols in each period's portfolio
				assets <- select(rets, symbol, weight) %>% 
					arrange(symbol)
				
				# Format returns tibble
				rets_df <- rets_tibble %>% 
					filter(symbol %in% assets$symbol) %>% 
					pivot_wider(names_from = symbol, values_from = ret.i) %>% 
					select(date, order(colnames(.))) %>%
					ungroup() %>% 
					distinct(date, .keep_all = T) %>% 
					arrange(date) %>% 
					select(-date)
				
				rets <- rets %>% 
					add_column(p.rho = private$.fn_p_rho(rets_df)) %>% 
					add_column(p.var = private$.fn_p_var(rets_df, assets$weight))
				
				portfolio <- rets %>% 
					mutate(set = "SingleIndex") %>% 
					select(
						set, symbol, exp.i, beta, std.i, weight, exp.p, p.beta, p.rho, p.var
					)
				
				private$.aa_output <- rets
				private$.aa_portfolio <- portfolio  
			} # end of asset_allocation()
		), # end of public list
		
		active = list(
			backtest_portfolios = function() {
				private$.backtest_portfolios
			}, # end of backtest_portfolios
			
			backtest_output = function() {
				private$.backtest_output
			}, # end of backtest_output
			
			aa_output = function() {
				private$.aa_output
			}, # end of aa_output
			
			aa_portfolio = function() {
				private$.aa_portfolio
			} # end of aa_portfolio
		), # end of active list
		
		lock_objects = F
	) # end of SingleIndex R6 class

# Define the R6 class for the 'Constant Correlation' optimization approach
ConsCorr <- 
	R6Class(
		"ConsCorr",
		private = list(
			
			# A function to model betas
			.fn_si = function(df) {
				lm(ret.i ~ ret.m, data = df, na.action = "na.exclude")
			},
			
			# A function to calculate portfolio variance
			.fn_p_var = function(rets_df, weights) {
				rets_df %>% 
					cov(., use = "pairwise.complete.obs") %>% 
					as.vector() %>% 
					`*`(., rep(weights, self$n_pos)) %>% 
					`*`(., rep(weights, each = self$n_pos)) %>% 
					sum(., na.rm = T)
			},
			
			# A function to calculate portfolio correlation
			.fn_p_rho = function(rets_df) {
				rets_df %>% 
					cor(., use = "pairwise.complete.obs") %>% 
					.[lower.tri(.)] %>% 
					mean(., na.rm = T)
			},
			
			.backtest_portfolios = NA,
			.backtest_output = NA,
			.aa_portfolio = NA,
			.aa_output = NA
		), # end of private list
		
		public = list(
			initialize = function(
				r_tibble = NULL, n_pos = NULL, roll_win = NULL, ...
			) {
				self$r_tibble <- r_tibble
				self$n_pos <- n_pos
				self$roll_win <- roll_win
			}, # end of initialize()
			
			r_tibble = NULL,
			n_pos = 10000,
			roll_win = 52,
			
			backtest = function() {
				
				# Use the opportunity set's mean rho as the constant correlation 
				# coefficient
				rho <- self$r_tibble %>% 
					select(symbol, date, ret.i) %>% 
					pivot_wider(names_from = symbol, values_from = ret.i) %>% 
					arrange(date) %>%  
					select(-date) %>% 
					private$.fn_p_rho()
				
				# Set the rolling functions
				rolling_si <- rollify(.f = function(ret.i, ret.m) {
					lm(ret.i ~ ret.m)
				},
				window = self$roll_win,
				unlist = F
				)
				
				rolling_mean <- rollify(.f = ~mean(.x, na.rm = T),
																window = self$roll_win
				)
				
				rolling_sd <- rollify(.f = ~sd(.x, na.rm = T),
															window = self$roll_win
				)
				
				# Run the backtest
				rets <- self$r_tibble %>% 
					group_by(symbol) %>% 
					
					# For rolling functions, observations must be at least as big as the 
					# rolling window
					filter(n() >= self$roll_win) %>% 
					mutate(
						exp.i = rolling_mean(ret.i), # Asset expected return
						std.i = rolling_sd(ret.i), # Asset standard deviation of returns
						std.m = rolling_sd(ret.m), # Market standard deviation for comparis.
						roll.si = rolling_si(ret.i, ret.m) # To later get asset beta
					) %>% 
					filter(!is.na(exp.i)) %>% # Remove unnecessary rows
					mutate(tidied = map(roll.si, tidy)) %>%
					unnest(tidied) %>%
					filter(term == "ret.m")%>% # keep only the beta term
					rename(beta = estimate) %>% 
					select(
						symbol, date, ret.i, ret.m, ret.f, exp.i, std.i, beta, std.m
					) %>% 
					mutate(excess = (exp.i - ret.f) / std.i) %>% 
					.[!is.infinite(.$excess), ] %>% 
					arrange(desc(date)) %>% 
					ungroup() %>% 
					group_by(date) %>% 
					arrange(desc(excess), .by_group = T) %>% 
					
					# Compute the cutoff coefficient
					mutate(
						id = seq_len(n()),
						rho.coef = (rho / (1 - rho + id * rho)), 
						cum.ex = cumsum(excess),
						c.coef = rho.coef * cum.ex
					)
				
				# Extract the cutoff rates
				c.star <- rets %>% 
					filter(excess >= c.coef, .preserve = T) %>% 
					top_n(1, c.coef) %>% 
					rename(c.star = c.coef) %>% 
					select(date, c.star)
				
				# Apply the cutoff rate
				rets <- rets %>% 
					left_join(c.star, by = "date") %>% 
					mutate(c.star = replace_na(c.star, 0),
								 
								 # Compute the investment proportion coefficient
								 z.coef = ((1 - rho) * std.i)^-1 * (excess - c.star)
					) %>% 
					
					# Allowing for short sales, rank assets by profitability potential
					arrange(desc(abs(z.coef)), .by_group = T) %>% 
					
					# Apply the cardinality constraint
					top_n(self$n_pos, abs(z.coef)) %>% 
					
					# Apply Lintner's short sale definition
					mutate(
						weight = z.coef / sum(abs(z.coef)),
						p.beta = sum(abs(weight) * beta)
					) %>% 
					arrange(symbol, .by_group = T) %>% 
					ungroup()
				
				# Compute portfolio variances and mean correlations
				# Get the symbols in each period's portfolio
				asset_names <- rets %>% 
					select(date, symbol) %>% 
					pivot_wider(names_from = symbol, values_from = symbol) %>% 
					ungroup() %>% 
					arrange(date) %>% 
					select(date, order(colnames(.)))
				
				# Format returns tibble
				rets_df <- self$r_tibble %>% 
					select(symbol, date, ret.i) %>% 
					pivot_wider(names_from = symbol, values_from = ret.i) %>% 
					select(date, order(colnames(.))) %>%
					ungroup() %>% 
					distinct(date, .keep_all = T) %>% 
					arrange(date)
				
				# Generate the rolling windows
				idx <- map(
					seq(1, nrow(rets_df) - self$roll_win, by = 1), 
					function(i) seq(i, i + self$roll_win - 1)
				)
				
				# Get the returns data frame that reflects each period's optimized 
				# portfolio
				roll_rets <- map(idx, function(i) rets_df[i, ])
				
				roll_rets <- map(
					seq_along(roll_rets),
					function(i)
						roll_rets[[i]][, colnames(roll_rets[[i]]) %in% asset_names[i, ]]
				)
				
				# Generate the variance-covariance matrices
				roll_cov <- map(
					seq_along(roll_rets),
					function(i)
						cov(roll_rets[[i]], use = "pairwise.complete.obs")
				)
				
				# Get corresponding weights
				asset_weights <- rets$weight %>% 
					matrix(., ncol = self$n_pos, byrow = T)
				
				# Calculate portfolio variance
				p_var <- map_dbl(
					seq_along(roll_cov),
					function(i) {
						as.vector(roll_cov[[i]]) %>% 
							`*`(., rep(asset_weights[i, ], self$n_pos)) %>% 
							`*`(., rep(asset_weights[i, ], each = self$n_pos)) %>% 
							sum(., na.rm = T)
					}
				)
				
				# Generate the correlation matrices
				roll_rho <- map(
					seq_along(roll_rets),
					function(i)
						cor(roll_rets[[i]], use = "pairwise.complete.obs")
				)
				
				# Calculate mean portfolio correlations
				p_rho <- map_dbl(
					seq_along(roll_rho),
					function(i) {
						roll_rho[[i]][lower.tri(roll_rho[[i]])] %>% mean(., na.rm = T)
					}
				)
				
				rets <- rets %>% 
					add_column(
						p.var = c(
							rep(p_var, each = self$n_pos),
							if(length(p_var) * self$n_pos < nrow(rets))
								rep(NA, self$n_pos)
						)
					) %>% 
					add_column(
						p.rho = c(
							rep(p_rho, each = self$n_pos),
							if(length(p_rho) * self$n_pos < nrow(rets))
								rep(NA, self$n_pos)
						)
					)
				
				# Compute portfolio performance and format appropriately
				portfolios <- rets %>% 
					group_by(symbol) %>% 
					
					# A rolling window backtest performance is computed with lagged 
					# weights
					mutate(
						ret.b = lag(weight, order_by = date) * ret.i,
						ret.b = replace_na(ret.b, 0), # Helpful for graphing purposes
						m.var = std.m^2
					) %>% 
					ungroup() %>% 
					group_by(date) %>% 
					mutate(ret.p = sum(ret.b)) %>% 
					
					# Format the portfolio tibble
					arrange(date) %>% 
					distinct(date, .keep_all = T) %>% 
					select(date, ret.p, p.beta, p.var, p.rho, ret.m, m.var) %>% 
					rename(ConsCorr = ret.p, index = ret.m) %>% 
					gather(index, ConsCorr, key = "set", value = "return") %>% 
					ungroup() %>% 
					group_by(set) %>% 
					
					# Compute compound returns
					mutate(
						return = replace(return, row_number() == 1, 0),
						level = cumprod(1 + return),
						p.beta = lag(p.beta),
						p.var = lag(p.var),
						p.rho = lag(p.rho)
					) %>% 
					ungroup()
				
				private$.backtest_portfolios  <- portfolios
				private$.backtest_output <- rets
			}, # end of backtest()
			
			asset_allocation = function() {
				
				# Get the last period risk-free rate
				rf_last <- self$r_tibble %>% 
					filter(date == max(date)) %>% 
					distinct(ret.f) %>% 
					as.numeric(.)
				
				# Run the asset allocation
				rets_tibble <- self$r_tibble %>% 
					select(symbol, date, ret.i) %>% 
					pivot_wider(names_from = symbol, values_from = ret.i) %>% 
					arrange(desc(date)) %>% 
					top_n(self$roll_win, wt = date) %>% 
					pivot_longer(2:ncol(.), names_to = "symbol", values_to = "ret.i") %>% 
					drop_na() %>% 
					group_by(symbol) %>% 
					arrange(date, .by_group = T) %>% 
					ungroup()
				
				to_join <- self$r_tibble %>% 
					select(date, ret.m, ret.f) %>% 
					distinct(date, .keep_all = T)
				
				rets <- rets_tibble %>% 
					group_by(symbol) %>% 
					left_join(to_join, by = "date") %>% 
					filter(n() > 1) %>% 
					ungroup()
				
				# Use the opportunity set's mean rho as the constant correlation 
				# coefficient
				rho <- rets %>% 
					select(symbol, date, ret.i) %>% 
					pivot_wider(names_from = symbol, values_from = ret.i) %>% 
					arrange(date) %>%  
					select(-date) %>% 
					private$.fn_p_rho()
				
				betas <- rets %>% 
					group_by(symbol) %>% 
					do(tidy(lm(ret.i ~ ret.m, .))) %>% 
					filter(term == "ret.m") %>% 
					rename(beta = estimate) %>% 
					select(symbol, beta) %>% 
					ungroup()
				
				rets <- rets %>% 
					add_column(rho = rho) %>% 
					left_join(betas, by = "symbol") %>% 
					group_by(symbol) %>% 
					mutate(
						exp.i = mean(ret.i, na.rm = T), # Asset expected return
						std.i = sd(ret.i, na.rm = T), # Asset standard deviation of returns
						ret.f = rf_last
					) %>% 
					filter(!is.na(exp.i)) %>% # Remove unnecessary rows
					select(
						symbol, ret.i, ret.m, ret.f, exp.i, std.i, beta
					) %>% 
					distinct(symbol, .keep_all = T) %>% 
					mutate(excess = (exp.i - ret.f) / std.i) %>% 
					.[!is.infinite(.$excess), ] %>% 
					ungroup() %>%  
					arrange(desc(excess)) %>%
					
					# Compute the cutoff coefficient
					mutate(
						id = seq_len(n()),
						rho.coef = (rho / (1 - rho + id * rho)), 
						cum.ex = cumsum(excess),
						c.coef = rho.coef * cum.ex
					)
				
				# Extract the cutoff rates
				c.star <- rets %>% 
					filter(excess >= c.coef, .preserve = T) %>% 
					top_n(1, c.coef) %>% 
					rename(c.star = c.coef) %>% 
					select(c.star)
				
				# Apply the cutoff rate
				rets <- rets %>% 
					add_column(c.star = c.star$c.star) %>% 
					
					# Compute the investment proportion coefficient
					mutate(z.coef = ((1 - rho) * std.i)^-1 * (excess - c.star)) %>% 
					
					# Allowing for short sales, rank assets by profitability potential
					arrange(desc(abs(z.coef))) %>% 
					
					# Apply the cardinality constraint
					top_n(self$n_pos, abs(z.coef)) %>% 
					
					# Apply Lintner's short sale definition
					mutate(
						weight = z.coef / sum(abs(z.coef)),
						exp.p = sum(weight * exp.i),
						p.beta = sum(abs(weight) * beta)
					) %>% 
					arrange(symbol) %>% 
					ungroup()
				
				# Compute portfolio variances and mean correlations
				# Get the symbols in each period's portfolio
				assets <- select(rets, symbol, weight) %>% 
					arrange(symbol)
				
				# Format returns tibble
				rets_df <- rets_tibble %>% 
					select(symbol, date, ret.i) %>% 
					filter(symbol %in% assets$symbol) %>% 
					pivot_wider(names_from = symbol, values_from = ret.i) %>% 
					select(date, order(colnames(.))) %>%
					ungroup() %>% 
					distinct(date, .keep_all = T) %>% 
					arrange(date) %>% 
					select(-date)
				
				rets <- rets %>% 
					add_column(p.rho = private$.fn_p_rho(rets_df)) %>% 
					add_column(p.var = private$.fn_p_var(rets_df, assets$weight))
				
				portfolio <- rets %>% 
					mutate(set = "ConsCorr") %>% 
					select(
						set, symbol, exp.i, beta, std.i, weight, exp.p, p.beta, p.rho, p.var
					)
				
				private$.aa_output <- rets
				private$.aa_portfolio <- portfolio  
			} # end of asset_allocation()
		), # end of public list
		
		active = list(
			backtest_portfolios = function() {
				private$.backtest_portfolios
			}, # end of backtest_portfolios
			
			backtest_output = function() {
				private$.backtest_output
			}, # end of backtest_output
			
			aa_portfolio = function() {
				private$.aa_portfolio
			}, # end of aa_portfolio
			
			aa_output = function() {
				private$.aa_output
			} # end of aa_output
		), # end of active list
		
		lock_objects = F
	) # end of ConsCorr R6 class

# Define the R6 class for the 'Minimum Variance with Target Return' 
# optimization approach
MinVar <- 
	R6Class(
		"MinVar",
		private = list(
			.p_0 = list(),
			.exp_r = NULL,
			.rets = NULL,
			.varcov = NULL,
			.rho = NULL,
			.curr_per = NULL,
			.m_var = NULL,
			
			# A function to estimate betas
			.fn_si = function(df) {
				lm(ret.i ~ ret.m, data = df, na.action = "na.exclude")
			},
			
			# A function to calculate portfolio variance
			.fn_p_var = function(rets_df, weights) {
				rets_df %>% 
					cov(., use = "pairwise.complete.obs") %>% 
					as.vector() %>% 
					`*`(., rep(weights, self$n_pos)) %>% 
					`*`(., rep(weights, each = self$n_pos)) %>% 
					sum(., na.rm = T)
			},
			
			# A function to calculate portfolio correlation
			.fn_p_rho = function(rets_df) {
				rets_df %>% 
					cor(., use = "pairwise.complete.obs") %>% 
					.[lower.tri(.)] %>% 
					mean(., na.rm = T)
			},
			
			.backtest_portfolios = NA,
			.backtest_output = NA,
			.aa_portfolio = NA,
			.aa_output = NA
		),
		public = list(
			initialize = function(
				r_tibble = NULL, n_pos = NULL, roll_win = NULL, target_ret = NULL, ...
			) {
				self$r_tibble <- r_tibble
				self$n_pos <- n_pos
				self$roll_win <- roll_win
				self$target_ret <- target_ret
			}, # end of initialize()
			
			r_tibble = NULL,
			n_pos = 10000,
			roll_win = 52,
			n_gen = 100,
			trar_p = 0.5,
			mut_p = 0.01,
			target_ret = 0.1/52,
			
			backtest = function() {
				
				# Turn the tibble data into a manageable xts object
				xts_rets <- self$r_tibble %>% 
					select(symbol, date, ret.i) %>% 
					pivot_wider(names_from = symbol, values_from = ret.i) %>%
					select(date, order(colnames(.))) %>% 
					arrange(date) %>% 
					tk_xts()
				
				# Create a dates list
				dates <- self$r_tibble %>%
					distinct(date) %>%
					arrange(date) %>% 
					top_n(nrow(.) - self$roll_win + 1) %>%
					as.matrix()
				
				# Generate the rolling windows
				idx <- lapply(
					seq(1, nrow(xts_rets) - self$roll_win, by = 1), 
					function(i) seq(i, i + self$roll_win - 1)
				)
				
				# Create the returns data frames for each period
				roll_rets <- lapply(idx, function(i) xts_rets[i, ])
				roll_rets <- lapply(
					seq_along(roll_rets),
					function(i)
						roll_rets[[i]][, colSums(!is.na(roll_rets[[i]])) >= 2]
				)
				
				# Generate the rolling variance-covariance and correlation matrices, as
				# well as the rolling expected returns
				roll_cov <- lapply(
					seq_along(roll_rets),
					function(i) {
						cov(roll_rets[[i]], use = "pairwise.complete.obs")
					}
				)
				
				roll_rho <- lapply(
					seq_along(roll_rets),
					function(i)
						cor(roll_rets[[i]], use = "pairwise.complete.obs")
				)
				
				roll_exp <- lapply(
					seq_along(roll_rets),
					function(i)
						colMeans(roll_rets[[i]], na.rm = T)
				)
				
				# Initialize backtest output tibble
				private$.backtest_output <- tibble(
					date = as_date(dates[1]),
					set = NA,
					symbol = NA,
					weight = NA,
					p.var = NA,
					exp.p = 0,
					ret.m = 0,
					beta = NA,
					p.rho = NA,
					m.var = NA
				)
				
				# Inputs for quadratic optimization (quadprog::solve.QP())
				bvec <-  c(1, self$target_ret)
				dvec <-  rep(0, self$n_pos)
				meq <-  1
				
				# Wrap the genetic algorithm in a function
				fn_ga <- function(t) {
					# Get this period's relevant data
					varcov <- roll_cov[[t]] 
					to_keep <- colnames(varcov)[colSums(is.na(varcov)) == 0]
					
					private$.varcov <- 
						varcov[row.names(varcov) %in% to_keep, colnames(varcov) %in% to_keep]
					
					rets <- roll_rets[[t]]
					private$.rets <-rets[, colnames(rets) %in% to_keep]
					
					private$.m_var <- var(roll_rets[[1]]$GSPC, na.rm = T) %>% as.numeric()
					
					rho <- roll_rho[[t]]
					private$.rho <- 
						rho[row.names(rho) %in% to_keep, colnames(rho) %in% to_keep]
					
					exp_r <- roll_exp[[t]][names(roll_exp[[t]]) %in% to_keep]
					private$.curr_per <- as_date(dates[t + 1])
					
					cat(paste(private$.curr_per, "computations"), "\n")
					
					# Initialize p_0 list of initial population
					p_0 <- list()
					p_0 <- lapply(
						seq_len(self$n_gen),
						function(i) {
							mu <- sample(exp_r, self$n_pos)
							mu <- mu[order(names(mu))]
							
							# Ensure variance-covariance matrix is positive definite for the
							# quadratic solver - revisit this solution
							varcov <- nearPD(
								private$.varcov[row.names(private$.varcov) %in% names(mu), names(mu)]
							)
							
							Amat <- cbind(rep(1, self$n_pos), mu)
							
							# solv.QP optimizes 1/2*w*Q*w^T, so to get the portfolio variance, 
							# the covariance matrix is multiplied by 2
							x <- solve.QP(2 * varcov$mat, dvec, Amat, bvec, meq)
							
							names(x$solution) <- names(mu)
							p_0[[i]] <- x
						}
					)
					
					p_0 <- Filter(Negate(anyNA), p_0)
					
					p_0 <- map(seq_along(p_0), function(p) {
						enframe(p_0[[p]]$solution) %>% 
							mutate(
								name = names(p_0[[p]]$solution),
								p.var = p_0[[p]]$value
							) %>% 
							add_column(set = paste0("p", p), .before = "name") %>% 
							rename(symbol = name, weight = value)
					}
					) %>% 
						enframe() %>% 
						unnest(value) %>% 
						select(-name) %>% 
						group_by(set) %>% 
						mutate(weight = weight / sum(abs(weight))) %>% 
						ungroup()
					
					private$.exp_r <- exp_r %>% 
						enframe() %>% 
						rename(symbol = name, exp.i = value)
					
					private$.p_0 <- left_join(p_0, private$.exp_r, by = "symbol")
					
					# Go to fn_transrar() from here
					fn_transrar()
				}
				
				# The genetic algorithm TransRAR crossover function
				fn_transrar <- function() { 
					
					# Run the TransRAR crossover for the current period
					map(seq_len(self$n_gen), function(x) {
						cat(paste(private$.curr_per, "generation", x), "\n")
						
						# Set the multiset-union
						parents <- private$.p_0 %>% 
							group_by(set) %>% 
							nest() %>% 
							ungroup() %>% 
							sample_n(2, replace = F) %>% 
							unnest(data) %>% 
							select(-p.var) %>% 
							spread(set, weight) %>% 
							rename(eye.1 = 2, eye.2 = 3) %>% 
							mutate(
								trar.p = case_when(
									!is.na(eye.1 & eye.2) ~ 1,
									!is.na(eye.1 | eye.2) ~ runif(nrow(.))
								)
							) %>% 
							select(symbol, trar.p) %>% 
							filter(trar.p >= self$trar_p)
						
						# Create the child chromosome
						child <- parents %>%
							sample_n(min(self$n_pos, nrow(parents)))
						
						# Remaining gene pool
						remainder <- private$.exp_r %>% 
							filter(!symbol %in% child$symbol) 
						
						# Fill child chromosome, if not complete
						gene_fill <- remainder %>% 
							sample_n(self$n_pos - nrow(child), replace = F)
						
						child <- child %>% 
							bind_rows(gene_fill) 
						
						# Apply mutation
						mutation <- runif(1)
						if(mutation < self$mut_p) {
							child <- child %>%
								sample_n(self$n_pos - 1) %>%
								bind_rows(sample_n(remainder, 1))
						}
						
						child <- child %>% 
							add_column(
								set = paste0("p", self$n_gen + x),
								.before = "symbol"
							) %>% 
							left_join(private$.exp_r, by = "symbol") %>% 
							rename(exp.i = exp.i.y) %>% 
							select(set, symbol, exp.i) %>% 
							arrange(symbol)
						
						# Calculate the child's fitness
						# Ensure variance-covariance matrix is positive definite for the
						# quadratic solver - revisit this solution
						varcov <- 
							nearPD(
								private$.varcov[row.names(private$.varcov) %in% child$symbol, child$symbol]
							)
						
						Amat <- cbind(rep(1, self$n_pos), child$exp.i)
						
						# solv.QP optimizes 1/2*w*Q*w^T, so to get the portfolio variance,  
						# the covariance matrix is multiplied by 2
						child_fit <- solve.QP(2 * varcov$mat, dvec, Amat, bvec, meq)
						child <- child %>% 
							mutate(
								weight = child_fit$solution / sum(abs(child_fit$solution)),
								p.var = child_fit$value
							) %>% 
							select(set, symbol, weight, p.var, exp.i)
						
						# Advance the n_gen portfolios with the least fitness value
						max_p_0_fit <- max(private$.p_0$p.var)
						if(child$p.var[1] < max_p_0_fit) {
							private$.p_0 <- private$.p_0 %>% 
								filter(p.var < max_p_0_fit) %>% 
								bind_rows(child)
						} 
					})
					
					private$.p_0 <- private$.p_0 %>% 
						group_by(set) %>% 
						mutate(
							weight = weight / sum(abs(weight)),
							exp.p = sum(exp.i * weight)
						) %>% 
						ungroup()
					
					# Choose portfolio for this period
					exp_m <- self$r_tibble %>% 
						filter(date == private$.curr_per) %>% 
						top_n(1, symbol) %>% 
						.$ret.m
					
					p_1 <- if (max(private$.p_0$exp.p < self$target_ret)) {
						
						# Consider adding a warning message to this regarding the out-of-
						# bounds target
						private$.p_0 %>% 
							filter(exp.p == max(exp.p)) %>% 
							add_column(date = private$.curr_per, .before = "set") %>% 
							mutate(
								ret.m = exp_m,
								m.var = private$.m_var
							)
					} else {
						private$.p_0 %>% 
							filter(exp.p >= self$target_ret) %>% 
							filter(p.var == min(p.var)) %>% 
							add_column(date = private$.curr_per, .before = "set") %>% 
							mutate(
								ret.m = exp_m,
								m.var = private$.m_var
							)
					}
					
					# Compute betas
					betas <- private$.rets %>% 
						as_tibble() %>% 
						add_column(ret.m = .$GSPC) %>% 
						select(p_1$symbol, ret.m) %>% 
						pivot_longer(
							1:(ncol(.) - 1), names_to = "symbol", values_to = "ret.i"
						) %>% 
						arrange(symbol) %>% 
						group_by(symbol) %>% 
						nest() %>% 
						mutate(model = map(data, private$.fn_si)) %>% 
						mutate(tidied = map(model, tidy)) %>% 
						unnest(tidied) %>% 
						filter(term == "ret.m") %>%  
						select(symbol, estimate) %>% 
						rename(beta = estimate) %>% 
						ungroup()
					
					# Compute mean correlation
					p.rho <- 
						private$.rho[row.names(private$.rho) %in% p_1$symbol, p_1$symbol]
					
					p.rho <- p.rho[lower.tri(p.rho)] %>% mean(., na.rm = T)
					
					p_1 <- p_1 %>% 
						left_join(betas, by = "symbol") %>% 
						add_column(p.rho = p.rho)
					
					# Add output to backtest tibble
					private$.backtest_output <- private$.backtest_output %>% 
						bind_rows(p_1)
					
					cat(paste(private$.curr_per, "computations completed"), "\n")
				} # end of fn_transrar
				
				map(seq_along(idx), fn_ga)
				
				private$.backtest_output <- private$.backtest_output %>% 
					left_join(
						self$r_tibble %>% 
							select(symbol, date, ret.i),
						by = c("symbol" = "symbol", "date" = "date")
					) %>% 
					ungroup() %>% 
					group_by(date) %>% 
					mutate(
						return = sum(ret.i * weight, na.rm = T),
						p.beta = mean(beta),
						set = "MinVar"
					) %>% 
					ungroup()  
				
				# Apply the format uniform across optimization approaches
				p_b <- private$.backtest_output %>%  
					select(date, p.beta, p.var, p.rho, set, return, m.var) %>% 
					distinct(date, .keep_all = T) %>% 
					bind_rows(
						private$.backtest_output %>% 
							distinct(date, .keep_all = T) %>% 
							mutate(set = "index") %>% 
							select(
								date, beta, p.var, p.rho, set, ret.m, m.var
							) %>% 
							rename(
								p.beta = beta, return = ret.m
							)
					) %>% 
					group_by(set) %>% 
					mutate(
						return = replace(return, row_number() == 1, 0),
						level = cumprod(1 + return)
					) %>% 
					ungroup()
				
				private$.p_0 <- list()
				private$.exp_r <- NULL
				private$.rets <- NULL
				private$.varcov <- NULL
				private$.rho <- NULL
				private$.curr_per <- NULL
				private$.m_var <- NULL
				
				private$.backtest_portfolios  <- p_b
			}, # end of backtest()
			
			asset_allocation = function() {
				
				# Run the asset allocation
				rets_tibble <- self$r_tibble %>% 
					select(symbol, date, ret.i) %>% 
					pivot_wider(names_from = symbol, values_from = ret.i) %>% 
					arrange(desc(date)) %>% 
					top_n(self$roll_win, wt = date) %>% 
					pivot_longer(2:ncol(.), names_to = "symbol", values_to = "ret.i") %>% 
					drop_na() %>% 
					group_by(symbol) %>% 
					arrange(date, .by_group = T) %>% 
					ungroup()
				
				to_join <- self$r_tibble %>% 
					select(date, ret.m, ret.f) %>% 
					distinct(date, .keep_all = T)
				
				rets_tibble <- rets_tibble %>% 
					group_by(symbol) %>% 
					left_join(to_join, by = "date") %>% 
					filter(n() > 1) %>% 
					ungroup()
				
				# Inputs for quadratic optimization (quadprog::solve.QP())
				bvec <-  c(1, self$target_ret)
				dvec <-  rep(0, self$n_pos)
				meq <-  1
				
				# The genetic algorithm TransRAR crossover function
				fn_transrar <- function() { 
					map(seq_len(self$n_gen), function(x) {
						cat(paste("generation", x), "\n")
						
						# Set the multiset-union
						parents <- private$.p_0 %>% 
							group_by(set) %>% 
							nest() %>% 
							ungroup() %>% 
							sample_n(2, replace = F) %>% 
							unnest(data) %>% 
							select(-c(p.var)) %>% 
							spread(set, weight) %>% 
							rename(eye.1 = 2, eye.2 = 3) %>% 
							mutate(
								trar.p = case_when(
									!is.na(eye.1 & eye.2) ~ 1,
									!is.na(eye.1 | eye.2) ~ runif(nrow(.))
								)
							) %>% 
							select(symbol, trar.p) %>% 
							filter(trar.p >= self$trar_p)
						
						# Create the child chromosome
						child <- parents %>%
							sample_n(min(self$n_pos, nrow(parents)))
						
						# Remaining gene pool
						remainder <- private$.exp_r %>% 
							filter(!symbol %in% child$symbol) 
						
						# Fill child chromosome, if not complete
						gene_fill <- remainder %>% 
							sample_n(self$n_pos - nrow(child), replace = F)
						
						child <- child %>% 
							bind_rows(gene_fill) 
						
						# Apply mutation
						mutation <- runif(1)
						if(mutation < self$mut_p) {
							child <- child %>%
								sample_n(self$n_pos - 1) %>%
								bind_rows(sample_n(remainder, 1))
						}
						
						child <- child %>% 
							add_column(
								set = paste0("p", self$n_gen + x),
								.before = "symbol"
							) %>% 
							left_join(private$.exp_r, by = "symbol") %>% 
							rename(exp.i = exp.i.y) %>% 
							select(set, symbol, exp.i) %>% 
							arrange(symbol)
						
						# Calculate the child's fitness
						# Ensure variance-covariance matrix is positive definite for the
						# quadratic solver - revisit this solution
						varcov <- 
							nearPD(
								private$.varcov[row.names(private$.varcov) %in% child$symbol, child$symbol]
							)
						
						Amat <- cbind(rep(1, self$n_pos), child$exp.i)
						child_fit <- solve.QP(2 * varcov$mat, dvec, Amat, bvec, meq)
						child <- child %>% 
							mutate(
								weight = child_fit$solution / sum(abs(child_fit$solution)),
								p.var = child_fit$value
							) %>% 
							select(set, symbol, weight, p.var, exp.i)
						
						# Advance the n_gen portfolios with the least fitness value
						max_p_0_fit <- max(private$.p_0$p.var)
						if(child$p.var[1] < max_p_0_fit) {
							private$.p_0 <- private$.p_0 %>% 
								filter(p.var < max_p_0_fit) %>% 
								bind_rows(child)
						} 
					})
					
					private$.p_0 <- private$.p_0 %>% 
						group_by(set) %>% 
						mutate(
							weight = weight / sum(abs(weight)),
							exp.p = sum(exp.i * weight)
						) %>% 
						ungroup()
					
					# Choose portfolio for this period
					market <- rets_tibble %>% 
						transmute(
							exp.m = mean(ret.m, na.rm = T),
							m.var = var(ret.m, na.rm = T)
						) %>% 
						distinct(exp.m, .keep_all = T)
					
					std.i <- rets_tibble %>% 
						group_by(symbol) %>% 
						transmute(std.i = sd(ret.i)) %>% 
						distinct(symbol, .keep_all = T)
					
					
					p_1 <- if (max(private$.p_0$exp.p < self$target_ret)) {
						
						# Consider adding a warning message to this regarding the out-of-
						# bounds target
						private$.p_0 %>% 
							filter(exp.p == max(exp.p)) %>% 
							mutate(
								exp.m = market$exp.m,
								m.var = market$m.var
							)
					} else {
						private$.p_0 %>% 
							filter(exp.p >= self$target_ret) %>% 
							filter(p.var == min(p.var)) %>% 
							mutate(
								exp.m = market$exp.m,
								m.var = market$m.var
							)
					}
					
					# Compute betas
					betas <- rets %>% 
						as_tibble() %>% 
						add_column(ret.m = .$GSPC) %>% 
						select(p_1$symbol, ret.m) %>% 
						pivot_longer(
							1:(ncol(.) - 1), names_to = "symbol", values_to = "ret.i"
						) %>% 
						arrange(symbol) %>% 
						group_by(symbol) %>% 
						nest() %>% 
						mutate(model = map(data, private$.fn_si)) %>% 
						mutate(tidied = map(model, tidy)) %>% 
						unnest(tidied) %>% 
						filter(term == "ret.m") %>%  
						select(symbol, estimate) %>% 
						rename(beta = estimate) %>% 
						ungroup()
					
					# Compute mean correlation
					p_rho <- rets %>% 
						select(p_1$symbol) %>% 
						private$.fn_p_rho()
					
					p_1 <- p_1 %>% 
						left_join(betas, by = "symbol") %>% 
						left_join(std.i, by = "symbol") %>% 
						add_column(p.rho = p_rho)
					
					# Add output to backtest tibble
					private$.aa_output <- p_1
				} # end of fn_transrar
				
				varcov <- rets_tibble %>% 
					select(symbol, date, ret.i) %>% 
					pivot_wider(names_from = symbol, values_from = ret.i) %>% 
					arrange(date) %>% 
					select(-date) %>% 
					cov(., use = "pairwise.complete.obs")
				
				to_keep <- colnames(varcov)[colSums(is.na(varcov)) == 0]
				
				private$.varcov <- 
					varcov[row.names(varcov) %in% to_keep, colnames(varcov) %in% to_keep]
				
				rets <- rets_tibble %>% 
					select(symbol, date, ret.i) %>% 
					pivot_wider(names_from = symbol, values_from = ret.i) %>% 
					arrange(date) %>% 
					select(-date)
				
				private$.rets <-rets[, colnames(rets) %in% to_keep]
				
				private$.m_var = var(rets$GSPC, na.rm = T) %>% as.numeric()
				
				exp_r <- rets %>% 
					colMeans(., na.rm = T) %>% 
					.[names(.) %in% to_keep]
				
				p_0 <- list()
				p_0 <- lapply(
					seq_len(self$n_gen),
					function(i) {
						mu <- sample(exp_r, self$n_pos)
						mu <- mu[order(names(mu))]
						
						# Ensure variance-covariance matrix is positive definite for the
						# quadratic solver - revisit this solution
						varcov <- nearPD(
							private$.varcov[row.names(private$.varcov) %in% names(mu), names(mu)]
						)
						
						Amat <- cbind(rep(1, self$n_pos), mu)
						x <- solve.QP(2 * varcov$mat, dvec, Amat, bvec, meq)
						
						names(x$solution) <- names(mu)
						p_0[[i]] <- x
					}
				)
				
				p_0 <- Filter(Negate(anyNA), p_0)
				
				p_0 <- map(seq_along(p_0), function(p) {
					enframe(p_0[[p]]$solution) %>% 
						mutate(
							name = names(p_0[[p]]$solution),
							p.var = p_0[[p]]$value
						) %>% 
						add_column(set = paste0("p", p), .before = "name") %>% 
						rename(symbol = name, weight = value)
				}
				) %>% 
					enframe() %>% 
					unnest(value) %>% 
					select(-name) %>% 
					group_by(set) %>% 
					mutate(weight = weight / sum(abs(weight))) %>% 
					ungroup()
				
				private$.exp_r <- exp_r %>% 
					enframe() %>% 
					rename(symbol = name, exp.i = value)
				
				private$.p_0 <- left_join(p_0, private$.exp_r, by = "symbol")
				
				# Go to fn_transrar() from here
				fn_transrar()
				
				# Apply the format uniform across optimization approaches
				p_b <- private$.aa_output %>%
					mutate(
						p.beta = sum(abs(weight) * beta),
						set = "MinVar"
					) %>% 
					select(
						set, symbol, weight, exp.i, std.i, beta, exp.p, p.var, p.rho, p.beta
					) 
				
				private$.p_0 <- list()
				private$.exp_r <- NULL
				private$.rets <- NULL
				private$.varcov <- NULL
				private$.rho <- NULL
				private$.curr_per <- NULL
				private$.m_var <- NULL
				
				private$.aa_portfolio  <- p_b
				
			} # end of asset_allocation()
		), # end of public list
		
		active = list(
			backtest_portfolios = function() {
				private$.backtest_portfolios
			}, # end of backtest_portfolios
			
			backtest_output = function() {
				private$.backtest_output
			}, # end of backtest_output
			
			aa_portfolio = function() {
				private$.aa_portfolio
			}, # end of aa_portfolio
			
			aa_output = function() {
				private$.aa_output
			} # end of aa_output
		), # end of active list
		
		lock_objects = F
	) # end of MinVar R6 class

# Initiate objects for the three classes
n_SingleIndex <- SingleIndex$new()
n_ConsCorr <- ConsCorr$new()
n_MinVar <- MinVar$new()


# ##############################################################################
# DATA PREPARATION AND CARDINALITY CONSTRAINT FUNCTIONS
# ##############################################################################

fn_tibble <- function(freq_string, dates_vec = c(min_date, max_date)) {
	r_tibble <- data %>% 
		filter(freq == freq_string) %>% 
		filter(date >= dates_vec[1] & date <= dates_vec[2]) %>% 
		left_join(filter(., symbol == "GSPC"), by = "date") %>% 
		left_join(filter(., symbol.x == "DGS10"), by = "date") %>% 
		rename(
			symbol = symbol.x.x, ret.i = ret.x.x, ret.m = ret.y.x, ret.f = ret.x.y
		) %>% 
		select(symbol, date, ret.i, ret.m, ret.f) %>% 
		filter(symbol != "DGS10") %>% 
		
		# The adjustment to weekly data creates irregular periodicity, which the
		# following filter solves for the purposes of this app
		group_by(date) %>% 
		filter(n() >= n_assets/2) %>% 
		ungroup() %>% 
		
		# Get rid of outliers that mess everything up
		group_by(symbol) %>% 
		filter(
			ret.i < mean(ret.i) + 2 * sd(ret.i) | ret.i > mean(ret.i) - 2 * sd(ret.i)
		) %>% 
		ungroup()
} # end of fn_tibble()

# Cardinality constraint application
fn_card <- function(card) {
	case_when(
		card == 0 ~ as.double(nrow(data)),
		card != 0 ~ as.double(card)
	)
}

# ##############################################################################
# EXPLORATORY DATA ANALYSIS FUNCTIONS
# ##############################################################################

# A helper function to filter the data for EDA
fn_eda_filter <- function(eda_freq, eda_dates) {
	data %>% 
		filter(freq == eda_freq) %>% 
		filter(date >= eda_dates[1] & date <= eda_dates[2]) %>% 
		filter(symbol != "DGS10")
} # end of fn_eda_filter 

# A function to populate eda_stats
fn_eda_stats <- function(eda_freq, eda_dates) {
	data_es <- fn_eda_filter(eda_freq, eda_dates) %>% 
		group_by(symbol) %>% 
		
		# Compute preliminary statistics
		mutate(
			mu = mean(ret),
			std.d = sd(ret)
		) %>% 
		
		# Prepare tibble to calculate betas
		left_join(filter(., symbol == "GSPC"), by = "date") %>% 
		rename(
			symbol = symbol.x, ret = ret.x, mu = mu.x, std.d = std.d.x, mkt = ret.y
		) %>% 
		select(symbol, date, ret, mu, std.d, mkt) 
	
	# Compute betas
	betas <- data_es %>% 
		do(tidy(lm(ret ~  mkt, .))) %>% 
		filter(term == "mkt") %>% 
		rename(beta = estimate) %>% 
		select(symbol, beta, p.value)
	
	data_es <- left_join(data_es, betas, by = "symbol") %>% ungroup()
	
	# Calculate mean correlation
	rho <- data_es %>%
		select(symbol, date, ret) %>%
		spread(symbol, ret) %>%
		select(-date) %>%
		cor(., use = "pairwise.complete.obs")
	
	rho <- rho[lower.tri(rho)] %>% mean(., na.rm = T)
	
	# Create the preliminary statistics table
	data_es <- data_es %>%
		distinct(symbol, .keep_all = T) %>%
		select(symbol, mu, std.d, beta)
	
	index_stats <- data_es %>%
		filter(symbol == "GSPC") %>%
		select(-symbol) %>%
		rename(mean_return = mu, mean_std_dev = std.d, mean_beta = beta) %>%
		gather("statistic", "index")
	
	data_es <- data_es %>%
		transmute(
			mean_return = mean(mu),
			mean_std_dev = mean(std.d),
			mean_beta = mean(beta)
		) %>%
		distinct(mean_beta, .keep_all = T) %>%
		gather("statistic", "invst_space") %>%
		left_join(index_stats, by = "statistic") %>%
		add_row(
			statistic = "mean_corr",
			invst_space = rho,
			index = 1
		)
	
	data_es
} # end of fn_eda_stats()

# A function to provide the data for eda_invst_mu_sigma
fn_eda_invst_mu_sigma <- function(eda_freq, eda_dates) {
	data_ms <- fn_eda_filter(eda_freq, eda_dates) %>% 
		group_by(symbol) %>% 
		
		# Compute mean return and standard deviation
		mutate(
			mu = mean(ret),
			std.d = sd(ret)
		) %>% 
		select(symbol, mu, std.d) %>% 
		distinct(symbol, .keep_all = T)
	
	data_ms
} # end of fn_eda_invst_mu_sigma()

# A function to provide the data for eda_invst_mu_beta
fn_eda_invst_mu_beta <- function(eda_freq, eda_dates) {
	data_mb <- fn_eda_filter(eda_freq, eda_dates) %>% 
		group_by(symbol) %>% 
		mutate(mu = mean(ret)) %>% 
		
		# Prepare tibble to calculate betas
		left_join(filter(., symbol == "GSPC"), by = "date") %>% 
		rename(
			symbol = symbol.x, ret = ret.x, mu = mu.x, mkt = ret.y
		) %>% 
		select(symbol, date, ret, mkt, mu) 
	
	# Compute betas
	betas <- data_mb %>% 
		do(tidy(lm(ret ~  mkt, .))) %>% 
		filter(term == "mkt") %>% 
		rename(beta = estimate) %>% 
		select(symbol, beta, p.value)
	
	data_mb <- data_mb %>% 
		left_join(betas, by = "symbol") %>% 
		ungroup() %>% 
		distinct(symbol, .keep_all = T) %>% 
		select(symbol, mu, beta, p.value)
	
	data_mb
} # end of fn_eda_invst_mu_beta()

# A function to provide the data for eda_invst_pdf_mu
fn_eda_invst_pdf_mu <- function(eda_freq, eda_dates) {
	data_pm <- fn_eda_filter(eda_freq, eda_dates) %>% 
		group_by(symbol) %>% 
		mutate(mu = mean(ret)) %>% 
		select(symbol, mu) %>% 
		distinct(symbol, .keep_all = T)
	
	data_pm
} # end of fn_eda_invst_pdf_mu()

# A function to provide the data for eda_invst_pdf_sigma
fn_eda_invst_pdf_sigma <- function(eda_freq, eda_dates) {
	data_ps <- fn_eda_filter(eda_freq, eda_dates) %>% 
		group_by(symbol) %>% 
		mutate(std.d = sd(ret)) %>% 
		select(symbol, std.d) %>% 
		distinct(symbol, .keep_all = T)
	
	data_ps
} # end of fn_eda_invst_pdf_sigma()

# A function to provide the data for eda_index_time
fn_eda_index_time <- function(eda_freq, eda_dates) {
	data_it <- data2 %>% 
		filter(date >= eda_dates[1] & date <= eda_dates[2]) %>% 
		filter(symbol == "GSPC")
	
	if (eda_freq != "daily") {
		data_it <- data_it %>% 
			tq_transmute_(
				select = "obs",
				mutate_fun = paste0("to.", eda_freq),
				indexAt = "lastof"
			)
	}
	
	data_it
} # end of fn_eda_index_time()

# The following function will be used for both eda_index_ret_time and 
# eda_index_pdf_ret
fn_eda_index_rets <- function(eda_freq, eda_dates) {
	data_irt <- fn_eda_filter(eda_freq, eda_dates) %>% 
		filter(symbol == "GSPC")
	
	data_irt
} # end of fn_eda_index_rets()


# ##############################################################################
# BACKTEST FUNCTIONS
# ##############################################################################

# A function to populate b_stats
fn_b_stats <-  function(b_optim) {
	b_stats <- 
		eval(
			parse(text = paste0("n_", b_optim[1], "$backtest_portfolios"))
		) %>%
		filter(set == "index") %>% 
		
		# Calculate statistics for the table
		mutate(
			mean_ret = mean(return, na.rm = T),
			mean_std_dev = mean(sqrt(m.var), na.rm = T), 
			mean_beta = 1,
			mean_corr = 1
		) %>% 
		
		# Apply an appropriate format
		distinct(mean_ret, .keep_all = T) %>% 
		select(mean_ret, mean_std_dev, mean_beta, mean_corr) %>% 
		add_column(approach = "Index", .before = "mean_ret")
	
	b_stats <- b_stats %>% 
		bind_rows(
			map_dfr(
				b_optim,
				function(x) {
					eval(parse(text = paste0("n_", x, "$backtest_portfolios"))) %>%
						filter(set == x) %>% 
						
						# Calculate statistics for the table
						mutate(
							mean_ret = mean(return, na.rm = T),
							mean_std_dev = mean(sqrt(p.var), na.rm = T),
							mean_beta = mean(p.beta, na.rm = T),
							mean_corr = mean(p.rho, na.rm = T)
						) %>% 
						
						# Apply an appropriate format
						distinct(mean_ret, .keep_all = T) %>% 
						select(mean_ret, mean_std_dev, mean_beta, mean_corr) %>% 
						add_column(approach = x, .before = "mean_ret") 
				}
			)
		)
	
	b_stats
} # end of fn_b_stats

# A function to provide data for backtest plots
fn_b_performance <- function(b_optim) {
	
	b_performance <- map_dfr(
		b_optim,
		function(x) {
			eval(parse(text = paste0("n_", x, "$backtest_portfolios"))) %>%
				filter(set == x) 
		}
	)
	
	b_performance <- eval(
		parse(text = paste0("n_", b_optim[1], "$backtest_portfolios"))
	) %>%
		filter(set == "index") %>%
		bind_rows(b_performance) %>% 
		mutate(p.std = sqrt(p.var))
	
	b_performance
} # end of fn_b_performance


# ##############################################################################
# ASSET ALLOCATION
# ##############################################################################

fn_aa_performance <- function(aa_optim) {
	aa_performance <- map_dfr(
		aa_optim,
		function(x) {
			eval(parse(text = paste0("n_", x, "$aa_portfolio"))) %>%
				filter(set == x) %>%
				mutate(std.p = sqrt(p.var)) %>% 
				select(
					set, symbol, weight, exp.i, std.i , beta, exp.p, std.p, p.beta, p.rho
				)
		}
	)
	
	aa_performance
} # end of fn_aa_performance()

fn_efficient_frontier <- function(aa_freq, aa_window, aa_optim) {
	eff_df <- fn_tibble(aa_freq) %>% 
		select(symbol, date, ret.i) %>% 
		pivot_wider(names_from = symbol, values_from = ret.i) %>% 
		arrange(desc(date)) %>% 
		top_n(aa_window, wt = date) %>% 
		arrange(date) %>% 
		select(-date)
	
	eff_df <- eff_df %>% 
		.[, colSums(., na.rm = T) > 0] %>% 
		.[, sum(!is.na(.)) > 2]
	
	eff_sd <- apply(eff_df, 2, sd, na.rm = T)
	eff_rets <- colMeans(eff_df, na.rm = T)
	
	min_sd <- min(eff_sd)
	names(min_sd) <- names(which.min(eff_sd))
	max_ret <- max(eff_rets)
	names(max_ret) <- names(which.max(eff_rets))
	
	min_sd_rets <- eff_df %>% 
		select(names(min_sd)) %>% 
		drop_na() %>% 
		mutate(id = seq_len(nrow(.)))
	
	max_ret_rets <- eff_df %>% 
		select(names(max_ret)) %>% 
		drop_na() %>% 
		mutate(id = seq_len(nrow(.)))
	
	eff_cov <- left_join(min_sd_rets, max_ret_rets, by = "id") %>% 
		select(-id) %>% 
		cov(., use = "pairwise.complete.obs")
	
	eff_f <- tibble(
		set = "Efficient_Frontier",
		lambda = seq(-0.5, 1, 0.001)
	) %>% 
		mutate(
			mu = min_sd * (1 - lambda) + max_ret * lambda,
			sigma = 
				(1 - lambda)^2 * eff_cov[1, 1] + 
				lambda^2 * eff_cov[2, 2] +
				2 * (1 - lambda) * lambda * eff_cov[1, 2]
		)
	
	optims <- fn_aa_performance(aa_optim) %>% 
		mutate(lambda = NA) %>% 
		select(set, lambda, exp.p, std.p) %>% 
		distinct(set, .keep_all = T) %>% 
		rename(mu = exp.p, sigma = std.p) %>% 
		bind_rows(eff_f)
		
	optims
} # end of fn_efficient_frontier


# ##############################################################################
# USER INTERFACE
# ##############################################################################

header <- dashboardHeader(disable = T) # end of dashboardHeader

sidebar <- dashboardSidebar(
	# Website title
	h2(strong("Classic Portfolio Analysis")),
	
	sidebarMenu(
		
		# Exploratory Data Analysis tab
		menuItem(
			"Exploratory Data Analysis", tabName = "eda"
		), # end of eda menuItem
		
		# Backtest tab
		menuItem(
			"Backtest", tabName = "backtest"
		), # end of backtest menuItem
		
		# Asset Allocation tab
		menuItem(
			"Asset Allocation", tabName = "asset_allocation"
		) # end of asset_allocation menuItem
	), # end of sidebarMenu
	
	# Data download
	radioButtons(
		"data_format",
		"Download Raw Data",
		choices = list(
			"Long form" = "data2asis",
			"Data frame" = "data2spread"
		) # end of choices data_format
	), # end of data_format radioButtons
	
	downloadButton(
		"data_download", 
		"Download"
	) # end of data_download downloadButon
) # end of dashboardSidebar

body <-	dashboardBody(
	tabItems(
		
		# Exploratory Data Analysis body
		tabItem(
			tabName = "eda",
			column(
				width = 3,
				
				# Exploratory Data Analysis parameter box
				box(
					width = 12,
					title = h3("Exploratory Data Analysis"),
					
					# EDA date range input
					dateRangeInput(
						"eda_dates",
						"Date Range",
						start = max_date - years(1),
						end = max_date,
						min = min_date,
						max = max_date,
						separator = "-"
					), # end of eda_dates dateRangeInput
					
					# EDA returns frequency
					selectInput(
						"eda_freq",
						"Returns Frequency",
						choices = list(
							"Daily" = "daily",
							"Weekly" = "weekly",
							"Monthly" = "monthly",
							"Quarterly" = "quarterly",
							"Yearly" = "yearly"
						), # end of choices eda_freq
						selected = "weekly"
					), # end of eda_freq radioButtons
					
					# EDA delayed action button
					actionButton(
						"eda_action",
						"Generate Plots"
					), # end of eda_action actionButton
					
					hr(),
					
					# Display preliminary statistics
					h3("Preliminary Statistics"),
					tableOutput("eda_stats")
				) # end of EDA parameter box
			), # end of EDA column 1 - the EDA box column
			
			# Investment space plots
			column(
				width = 5,
				
				h3("Investment Space"),
				
				plotOutput(
					"eda_invst_mu_sigma",
					height = eda_plot_heights,
					click = "eda1_click",
					dblclick = "eda1_dblclick",
					brush = brushOpts(
						id = "eda1_brush",
						resetOnNew = T
					)
				),
				
				h5("Clicked points"),
				tableOutput("eda1_click_info"),
				
				plotOutput(
					"eda_invst_mu_beta",
					height = eda_plot_heights,
					click = "eda2_click",
					dblclick = "eda2_dblclick",
					brush = brushOpts(
						id = "eda2_brush",
						resetOnNew = T
					)
				),
				
				h5("Clicked points"),
				tableOutput("eda2_click_info"),
				
				plotOutput("eda_invst_pdf_mu", height = eda_plot_heights),
				plotOutput("eda_invst_pdf_sigma", height = eda_plot_heights)
			), # end of column 2 - the investment space plots column
			
			# Index plots
			column(
				width = 4,
				
				h3("Index"),
				plotOutput("eda_index_time", height = eda_plot_heights),
				plotOutput("eda_index_ret_time", height = eda_plot_heights),
				plotOutput("eda_index_pdf_ret", height = eda_plot_heights),
			) # end of column 3 - the index plots column
		), # end of eda tabItem
		
		# Backtest body
		tabItem(
			tabName = "backtest",
			column(
				width = 3,
				
				# Backtest parameter box
				box(
					width = 12,
					title = h3("Backtest"),
					
					helpText(
						"Backtests are run as rolling windows over the date range selected."
					), # end of helpText in Backtest parameter box
					
					# Backtest date range input
					dateRangeInput(
						"b_dates",
						"Date Range",
						start = max_date - years(2),
						end = max_date,
						min = min_date,
						max = max_date,
						separator = "-"
					), # end of b_dates dateRangeInput
					
					# Backtest returns frequency
					selectInput(
						"b_freq",
						"Returns Frequency",
						choices = list(
							"Daily" = "daily",
							"Weekly" = "weekly",
							"Monthly" = "monthly",
							"Quarterly" = "quarterly",
							"Yearly" = "yearly"
						), # end of choices b_freq
						selected = "weekly"
					), # end of b_freq selectInput
					
					# Rolling window size
					numericInput(
						"b_window",
						"Rolling Window Size",
						52
					), # end of b_window numericInput
					
					helpText(
						"For weekly returns, a 52 window size reflects a year."
					), # end of Backtest rolling window size helpText
					
					# Cardinality constraint
					numericInput(
						"b_card",
						"Cardinality Constraint",
						10
					), # end of b_card numericInput
					
					helpText(
						"For no cardinality constraint, input the number 0."
					), # end of Backtest cardinality constraint helpText
					
					# Choice of optimization approaches to backtest
					checkboxGroupInput(
						"b_optim",
						"Optimization Approach",
						choiceNames = optim_names,
						choiceValues = optim_values
					), # end of b_optim checkboxGroupInput
					
					numericInput(
						"b_target_ret",
						"Target Return (%)",
						0.2
					), # end of b_target_ret numericInput
					
					helpText(
						"Target return applies only to the minimum variance choice, and is
						applied per return frequency (10% may be a reasonable target for
						yearly returns frequency, 0.2% reasonable for weekly returns)."
					), # end of Backtest target return helpText
					
					# Backtest delayed reaction button
					actionButton(
						"b_action",
						"Run backtests"
					) # end of b_action actionButton
				) # end of Backtest parameter box
			), # end of Backtest column 1 - Backtest parameter box
			
			column(
				width = 9,
				
				h3("Backtest Results"),
				
				# First row shows the statistics table and options to save the run
				fluidRow(
					
					# Backtest Results statistics box
					box(
						width = 9,
						title = "Backtest Statistics",
						tableOutput("b_stats")
					), # end of Backtest Results statistics box
					
					# Options to bookmark the backtest choices and to save reports
					box(
						width = 3,
						title = "Save Current Run",
						bookmarkButton()
					)
				), # end of fluidRow 1 Backtest Results
				
				# Second row shows one graph of backtest performance and one plotly of
				# portfolio constituents per backtest period
				fluidRow(
					column(
						width = 6,
						plotOutput("b_performance", height = eda_plot_heights),
						plotOutput("b_mu_time", height = eda_plot_heights),
						plotOutput("b_sigma_time", height = eda_plot_heights),
						plotOutput("b_beta_time", height = eda_plot_heights),
						plotOutput("b_rho_time", height = eda_plot_heights)
					),
					
					column(
						width = 6,
						plotOutput("b_pdf_mu", height = eda_plot_heights),
						plotOutput("b_pdf_sigma", height = eda_plot_heights),
						plotOutput("b_pdf_beta", height = eda_plot_heights),
						plotOutput("b_pdf_rho", height = eda_plot_heights)
					)
				) # end of fluidRow 2 Backtest Results
			) # end of backtest column 2 - Backtest outputs
		), # end of backtest tabItem
		
		# Asset Allocation body
		tabItem(
			tabName = "asset_allocation",
			column(
				width = 3,
				
				# Asset Allocation parameter box
				box(
					width = 12,
					title = h3("Asset Allocation"),
					
					# Number of observations
					numericInput(
						"aa_window",
						"Number of Observations",
						52
					), # end of b_window numericInput
					
					helpText(
						"The number of observations define the sample size to estimate the 
						parameters for portfolio optimization."
					),
					
					# Asset Allocation returns frequency
					selectInput(
						"aa_freq",
						"Returns Frequency",
						choices = list(
							"Daily" = "daily",
							"Weekly" = "weekly",
							"Monthly" = "monthly",
							"Quarterly" = "quarterly",
							"Yearly" = "yearly"
						), # end of choices b_freq
						selected = "weekly"
					), # end of aa_freq selectInput
					
					# Cardinality constraint
					numericInput(
						"aa_card",
						"Cardinality Constraint",
						10
					), # end of aa_card numericInput
					
					helpText(
						"For no cardinality constraint, input the number 0."
					), # end of Asset Allocation cardinality constraint helpText
					
					# Asset Allocation optimization approaches
					checkboxGroupInput(
						"aa_optim",
						"Optimization Approach",
						choiceNames = optim_names,
						choiceValues = optim_values
					), # end of aa_optim checkboxGroupInput
					
					numericInput(
						"aa_target_ret",
						"Target Return (%)",
						0.2
					), # end of aa_target_ret numericInput
					
					helpText(
						"Target return applies only to the minimum variance choice, and is
						applied per return frequency (10% may be a reasonable target for
						yearly returns frequency, 0.2% reasonable for weekly returns)."
					), # end of Asset Allocation target return helpText
					
					# Backtest delayed reaction button
					actionButton(
						"aa_action",
						"Run Allocations"
					) # end of aa_action actionButton
				) # end of asset_allocation parameters box
			), # end of asset_allocation column 1 - Asset Allocation parameters
			
			column(
				width = 9,
				
				h3("Asset Allocation Outputs"),
				
				# First row 1 - save button
				fluidRow(
					box(
						width = 3,
						title = "Save Current Run",
						bookmarkButton()
					)
				), # end of fluidRow 1 Asset Allocation save button
				
				# Second row - two columns, first graph, second tables
				fluidRow(
					
					# Column for the efficient frontier and the 3D Asset plot
					column(
						width = 9,
						
						plotOutput("efficient_frontier", height = eda_plot_heights),
						
						box(
							width = 12,
							title = "Portfolio Constituents",
							tableOutput("asset_allocation")
						)
					) # end of column 1 fluidRow 2 Asset Allocation Outputs
				) # end of fluidRow 2 Asset Allocation Outputs
			) # asset_allocation column 2 - Outputs
		) # end of asset_allocation tabItem
	) # end of tabItems
) # end of dashboardBody

ui <- dashboardPage(header, sidebar, body, skin = "black")


# ##############################################################################
# SERVER FUNCTION
# ##############################################################################

server <- function(input, output, session) {
	
	### DATA DOWNLOAD
	data_download <- reactive({
		
		# Data format selection
		switch(
			input$data_format,
			"data2asis" = data2,
			"data2spread" = spread(data2, symbol, obs)
		) # end of data_format switch
	}) # end of data_download reactive
	
	# Data download button output
	output$data_download <- downloadHandler(
		filename = "assets_data.csv",
		content = function(file) {
			write_csv(data_download(), file)
		} # end of content function data_download
	) # end of data_download downloadHandler
	
	### EXPLORATORY DATA ANALYSIS TAB
	
	observeEvent(input$eda_action, {
		output$eda_stats <- renderTable({
			fn_eda_stats(input$eda_freq, input$eda_dates)
		}, digits = 5) # end of eda_stats renderTable
		
		ranges <- reactiveValues(x = NULL, y = NULL)
		
		output$eda_invst_mu_sigma <- renderPlot({
			fn_eda_invst_mu_sigma(input$eda_freq, input$eda_dates) %>% 
				ggplot(aes(std.d, mu)) +
				geom_point(aes(color = symbol)) +
				geom_vline(xintercept = 0) +
				geom_hline(yintercept = 0) +
				coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = F) +
				labs(
					title = "Returns to Standard Deviation (drag + double-click to zoom)",
					x = "Standard Deviation",
					y = "Returns"
				) +
				theme(legend.position = "none")
		}) # end of eda_invst_mu_sigma renderPlot
		
		observeEvent(input$eda1_dblclick, {
			brush <- input$eda1_brush
			if (!is.null(brush)) {
				ranges$x <- c(brush$xmin, brush$xmax)
				ranges$y <- c(brush$ymin, brush$ymax)
			} else {
				ranges$x <- NULL
				ranges$y <- NULL
			}
		}) # end of eda1_dblclick observeEvent
		
		output$eda1_click_info <- renderTable({
			nearPoints(
				fn_eda_invst_mu_sigma(input$eda_freq, input$eda_dates),
				input$eda1_click,
				addDist = T
			)
		},
		digits = 5) # end of eda1_click_info renderPlot
		
		ranges2 <- reactiveValues(x = NULL, y = NULL)
		
		output$eda_invst_mu_beta <- renderPlot({
			fn_eda_invst_mu_beta(input$eda_freq, input$eda_dates) %>% 
				ggplot(aes(beta, mu)) +
				geom_point(aes(color = symbol)) +
				geom_vline(xintercept = 0) +
				geom_hline(yintercept = 0) +
				coord_cartesian(xlim = ranges2$x, ylim = ranges2$y, expand = F) +
				labs(
					title = "Returns to Beta (drag + double-click to zoom)",
					x = "Beta",
					y = "Returns"
				) +
				theme(legend.position = "none") 
		}) # end of eda_invst_mu_beta renderPlot
		
		observeEvent(input$eda2_dblclick, {
			brush <- input$eda2_brush
			if (!is.null(brush)) {
				ranges2$x <- c(brush$xmin, brush$xmax)
				ranges2$y <- c(brush$ymin, brush$ymax)
			} else {
				ranges2$x <- NULL
				ranges2$y <- NULL
			}
		}) # end of eda1_dblclick observeEvent
		
		output$eda2_click_info <- renderTable({
			nearPoints(
				fn_eda_invst_mu_beta(input$eda_freq, input$eda_dates),
				input$eda2_click,
				addDist = T
			)
		},
		digits = 5) # end of eda2_click_info renderPlot
		
		output$eda_invst_pdf_mu <- renderPlot({
			fn_eda_invst_pdf_mu(input$eda_freq, input$eda_dates) %>% 
				ggplot(aes(mu)) +
				geom_density(color = "darkgreen") +
				geom_vline(xintercept = 0) +
				labs(
					title = "Density of Returns",
					x = "Returns"
				) +
				theme(
					axis.text.y = element_blank(),
					axis.ticks.y = element_blank(),
					axis.title.y = element_blank()
				)
		}) # end of eda_invst_pdf_mu renderPlot
		
		output$eda_invst_pdf_sigma <- renderPlot({
			fn_eda_invst_pdf_sigma(input$eda_freq, input$eda_dates) %>% 
				ggplot(aes(std.d)) +
				geom_density(color = "darkgreen") +
				geom_vline(xintercept = 0) +
				labs(
					title = "Density of Standard Deviations",
					x = "Standard Deviations"
				) +
				theme(
					axis.text.y = element_blank(),
					axis.ticks.y = element_blank(),
					axis.title.y = element_blank()
				)
		}) # end of eda_invst_pdf_sigma renderPlot
		
		output$eda_index_time <- renderPlot({
			fn_eda_index_time(input$eda_freq, input$eda_dates) %>% 
				ggplot(aes(date, obs)) +
				geom_line(color = "purple4") +
				labs(
					title = "Index Level (S&P 500)"
					) +
				theme(
					legend.position = "none",
					axis.title = element_blank()
				)
		}) # end of eda_index_time renderPlot
		
		output$eda_index_ret_time <- renderPlot({
			fn_eda_index_rets(input$eda_freq, input$eda_dates) %>% 
				ggplot(aes(date, ret)) +
				geom_line(color = "purple4") +
				geom_hline(yintercept = 0) +
				labs(
					title = "Index Returns over Time",
					y = "Returns"
				) +
				theme(
					legend.position = "none",
					axis.title.x = element_blank()
				)
		}) # end of eda_index_ret_time renderPlot
		
		output$eda_index_pdf_ret <- renderPlot({
			fn_eda_index_rets(input$eda_freq, input$eda_dates) %>% 
				ggplot(aes(ret)) +
				geom_density(color = "purple4") +
				geom_vline(xintercept = 0) +
				labs(
					title = "Density of Index Returns",
					x = "Index Returns"
				) +
				theme(
					legend.position = "none",
					axis.title.y = element_blank(),
					axis.text.y = element_blank(),
					axis.ticks.y = element_blank()
				)
		}) # end of eda_index_pdf_ret renderPlot
	}) # end of eda_action observeEvent
	
	### BACKTEST TAB
	observeEvent(input$b_action, {
		map(
			input$b_optim,
			function(x) {
				eval(
					parse(
						text = paste0(
							"n_", x, "$r_tibble <- fn_tibble(
							input$b_freq, dates_vec = input$b_dates
							)"
						)
					)
				)
				
				eval(parse(text = paste0("n_", x, "$n_pos <- fn_card(input$b_card)")))
				
				eval(parse(text = paste0("n_", x, "$roll_win <- input$b_window")))
				
				eval(
					parse(
						text = paste0("n_", x, "$target_ret <- input$b_target_ret / 100")
					)
				)
				
				eval(parse(text = paste0("n_", x, "$backtest()")))
			}
		)
		
		output$b_stats <- renderTable({fn_b_stats(input$b_optim)}, digits = 5)
		
		output$b_performance <- renderPlot({
			fn_b_performance(input$b_optim) %>% 
				ggplot(aes(date, level)) +
				geom_line(aes(color = set)) +
				labs(
					title = "Optimization Approach Performances"
				) +
				theme(
					axis.title.x = element_blank(),
					axis.title.y = element_blank()
				)
		}) # end of b_performance renderPlot
		
		output$b_mu_time <- renderPlot({
			fn_b_performance(input$b_optim) %>%
				ggplot(aes(date, return)) +
				geom_line(aes(color = set)) +
				labs(
					title = "Returns"
				) +
				theme(
					axis.title = element_blank()
				)
		}) # end of b_mu_time renderPlot
		
		output$b_sigma_time <- renderPlot({
			fn_b_performance(input$b_optim) %>%
				ggplot(aes(date, p.std)) +
				geom_line(aes(color = set)) +
				labs(
					title = "Changes in Standard Deviations",
					y = "Portfolio Standard Deviations"
				) +
				theme(
					axis.title.x = element_blank()
				)
		}) # end of b_sigma_time renderPlot
		
		output$b_beta_time <- renderPlot({
			fn_b_performance(input$b_optim) %>%
				ggplot(aes(date, p.beta)) +
				geom_line(aes(color = set)) +
				labs(
					title = "Changes in Beta",
					y = "Portfolio Beta"
				) +
				theme(
					axis.title.x = element_blank()
				)
		}) # end of b_beta_time renderPlot
		
		output$b_rho_time <- renderPlot({
			fn_b_performance(input$b_optim) %>%
				ggplot(aes(date, p.rho)) +
				geom_line(aes(color = set)) +
				labs(
					title = "Changes in Correlation",
					y = "Portfolio Correlation"
				) +
				theme(
					axis.title.x = element_blank()
				)
		}) # end of b_rho_time renderPlot
		
		output$b_pdf_mu <- renderPlot({
			fn_b_performance(input$b_optim) %>%
				ggplot(aes(return)) +
				geom_density(aes(color = set)) +
				labs(
					title = "Density of Returns",
					x = "Portfolio Returns"
				) +
				theme(
					axis.title.y = element_blank(),
					axis.text.y = element_blank(),
					axis.ticks.y = element_blank()
				)
		}) # end of b_pdf_mu renderPlot
		
		output$b_pdf_sigma <- renderPlot({
			fn_b_performance(input$b_optim) %>%
				ggplot(aes(p.std)) +
				geom_density(aes(color = set)) +
				labs(
					title = "Density of Portfolio Standard Deviations",
					x = "Standard Deviations"
				) +
				theme(
					axis.title.y = element_blank(),
					axis.text.y = element_blank(),
					axis.ticks.y = element_blank()
				)
		}) # end of b_pdf_sigma renderPlot
		
		output$b_pdf_beta <- renderPlot({
			fn_b_performance(input$b_optim) %>%
				ggplot(aes(p.beta)) +
				geom_density(aes(color = set)) +
				labs(
					title = "Density of Betas",
					x = "Portfolio Beta"
				) +
				theme(
					axis.title.y = element_blank(),
					axis.text.y = element_blank(),
					axis.ticks.y = element_blank()
				)
		}) # end of b_pdf_mu renderPlot
		
		output$b_pdf_rho <- renderPlot({
			fn_b_performance(input$b_optim) %>%
				ggplot(aes(p.rho)) +
				geom_density(aes(color = set)) +
				labs(
					title = "Density of Correlations",
					x = "Portfolio Correlation"
				) +
				theme(
					axis.title.y = element_blank(),
					axis.text.y = element_blank(),
					axis.ticks.y = element_blank()
				)
		}) # end of b_pdf_mu renderPlot
	}) # end of b_action observeEvent
	
	### ASSET ALLOCATION TAB
	observeEvent(input$aa_action, {
		map(
			input$aa_optim,
			function(x) {
				eval(
					parse(
						text = paste0(
							"n_", x, "$r_tibble <- fn_tibble(input$aa_freq)"
						)
					)
				)
				
				eval(parse(text = paste0("n_", x, "$n_pos <- fn_card(input$aa_card)")))
				
				eval(parse(text = paste0("n_", x, "$roll_win <- input$aa_window")))
				
				eval(
					parse(
						text = paste0("n_", x, "$target_ret <- input$aa_target_ret / 100")
					)
				)
				
				eval(parse(text = paste0("n_", x, "$asset_allocation()")))
			}
		)
		
		output$efficient_frontier <- renderPlot({
			 ggplot() +
				geom_path(
					data = fn_efficient_frontier(
						input$aa_freq, input$aa_window, input$aa_optim
					) %>%
						filter(set == "Efficient_frontier"),
					aes(sigma, mu),
					size = 0.5
				) +
				geom_point(
					data = fn_efficient_frontier(
						input$aa_freq, input$aa_window, input$aa_optim
					) %>%
						filter(set != "Efficient_frontier"),
					aes(sigma, mu, color = set),
					size = 2
				)
		}) # end of efficient_frontier renderPlot
		
		output$asset_allocation <- renderTable({
			fn_aa_performance(input$aa_optim)
		}, digits = 5)
	}) # end of aa_action observeEvent
} # end of server function

shinyApp(ui, server, enableBookmarking = "url")
