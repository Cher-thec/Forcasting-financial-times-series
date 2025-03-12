# Import des bibliothèques nécessaires
library(fGarch)
library(quantmod)
library(dygraphs)
library(moments)
library(urca)
library(forecast)
library(tseries)
library(xts)

# Importation des données
symbol_gold <- "GC=F"
symbol_euro <- "EURUSD=X"

start <- as.Date("2022-01-01")
end <- as.Date("2024-08-31")

gold_data <- getSymbols(symbol_gold, src = "yahoo", from = start, to = end, auto.assign = FALSE)
euro_data <- getSymbols(symbol_euro, src = "yahoo", from = start, to = end, auto.assign = FALSE)

# Suppression des données manquantes
gold_data <- na.omit(gold_data)
euro_data <- na.omit(euro_data)

par(mfrow = c(3, 2))

plot(euro_data[,1],main="Open EURO/USD",col="blue")
plot(euro_data[,2],main="High EURO/USD",col="blue")
plot(euro_data[,3],main="Low EURO/USD",col="blue")
plot(euro_data[,4],main="Close EURO/USD",col="blue")
plot(euro_data[,5],main="Volume EURO/USD",col="blue")
plot(euro_data[,6],main="Adjusted EURO/USD",col="blue")


plot(gold_data[,1],main="Open GOLD/USD",col="gold")
plot(gold_data[,2],main="High GOLD/USD",col="gold")
plot(gold_data[,3],main="Low GOLD/USD",col="gold")
plot(gold_data[,4],main="Close GOLD/USD",col="gold")
plot(gold_data[,5],main="Volume GOLD/USD",col="gold")
plot(gold_data[,6],main="Adjusted GOLD/USD",col="gold")

# Prix de clôture
gold_cl <- Cl(gold_data)
euro_cl <- Cl(euro_data)

# Plot the Data
par(mfrow=c(1,2))

dygraph(gold_cl, main = "Prix de clôture de l'Or (Gold/USD)") %>%
  dyAxis("x", label = "Date") %>%
  dyAxis("y", label = "Prix (USD)") %>%
  dyOptions(colors = "gold")

dygraph(euro_cl, main = "Prix de clôture Euro/USD") %>%
  dyAxis("x", label = "Date") %>%
  dyAxis("y", label = "Prix (USD)") %>%
  dyOptions(colors = "blue")

# Division en ensemble d'entraînement (80% des données) et ensemble de test (20% des données)
train_size <- 0.8
train_length_gold <- floor(train_size * length(gold_cl))
train_length_euro <- floor(train_size * length(euro_cl))

# Séparation des séries pour l'analyse (entraînement)
train_gold <- head(gold_cl, train_length_gold)
test_gold <- tail(gold_cl, length(gold_cl) - train_length_gold)

train_euro <- head(euro_cl, train_length_euro)
test_euro <- tail(euro_cl, length(euro_cl) - train_length_euro)

# Statistiques descriptives sur les ensembles d'entraînement
euro_stats <- c(
  mean = mean(train_euro, na.rm = TRUE),       
  median = median(train_euro, na.rm = TRUE),   
  variance = var(train_euro, na.rm = TRUE),    
  skewness = skewness(train_euro, na.rm = TRUE), 
  kurtosis = kurtosis(train_euro, na.rm = TRUE)
)

gold_stats <- c(
  mean = mean(train_gold, na.rm = TRUE),        
  median = median(train_gold, na.rm = TRUE),    
  variance = var(train_gold, na.rm = TRUE),     
  skewness = skewness(train_gold, na.rm = TRUE), 
  kurtosis = kurtosis(train_gold, na.rm = TRUE)
)

euro_stats
gold_stats

# Vérification de la stationnarité sur l'ensemble d'entraînement
adf.test(as.numeric(train_euro))
kpss_euro <- ur.kpss(as.numeric(train_euro))
summary(kpss_euro)

adf.test(as.numeric(train_gold))
kpss_gold <- ur.kpss(as.numeric(train_gold))
summary(kpss_gold)


# Différentiation pour rendre stationnaires
euro_diff <- na.omit(diff(train_euro))
gold_diff <- na.omit(diff(train_gold))

adf.test(as.numeric(euro_diff))
adf.test(as.numeric(gold_diff))

#Autoregressive and Moving Average Models

acf(euro_diff)
pacf(euro_diff)
acf(gold_diff)
pacf(gold_diff)

# Modélisation ARIMA sur l'ensemble d'entraînement
euro_arima <- auto.arima(euro_diff)
gold_arima <- auto.arima(gold_diff, seasonal = FALSE)

summary(euro_arima)
summary(gold_arima)

# Initialisation des variables
best_aic <- Inf  # Pour stocker le meilleur AIC
best_model <- NULL  # Pour stocker le meilleur modèle
best_p <- 0  # Meilleur paramètre AR
best_q <- 0  # Meilleur paramètre MA

# Boucle pour tester différentes valeurs de p et q (par exemple de 1 à 5)
for (p in 1:5) {
  for (q in 1:5) {
    # Essayer de créer un modèle ARIMA(p, 0, q)
    model <- tryCatch({
      arima_model <- Arima(euro_diff, order = c(p, 0, q))  # Modèle ARIMA(p, 0, q)
      arima_model
    }, error = function(e) NULL)  # Gérer les erreurs et éviter l'arrêt en cas de problème
    
    # Si le modèle est valide et n'est pas ARIMA(0,0,0), comparer les AIC
    if (!is.null(model) && !(p == 0 && q == 0)) {
      aic_value <- AIC(model)
      
      # Comparer les AIC et garder le meilleur modèle
      if (aic_value < best_aic) {
        best_aic <- aic_value
        best_model <- model
        best_p <- p
        best_q <- q
      }
    }
  }
}

# Résultats
cat("Meilleur modèle ARIMA(", best_p, ",0,", best_q, ") avec AIC =", best_aic, "\n")
summary(best_model)

#Fit AR and ARMA models to the train series

euro_model1<- Arima(euro_cl, order = c(2, 1, 3))
summary(euro_model1)
gold_model<- Arima(gold_cl, order = c(2, 1, 2))

euro_model3<- Arima(test_euro, order = c(2, 1, 3))
summary(euro_model1)
gold_model3<- Arima(test_gold, order = c(2, 1, 2))

#Residual Analysis

#tsdiag(euro_model1)

checkresiduals(euro_model1)
checkresiduals(gold_model)


#Heteroscedasticity Testing
euro_res<-residuals(euro_model1)
gold_res<-residuals(gold_model)

library(FinTS)
ArchTest(euro_res, lag=20)
ArchTest(gold_res, lags = 20)


#ARCH and GARCH Models
#Fitting ARCH Models

euro_arch <- garchFit(~ garch(1, 0), data = euro_res, trace = FALSE)
summary(euro_arch)

gold_arch <- garchFit(~ garch(1, 0), data = gold_res, trace = FALSE)
summary(gold_arch)

#GARCH(1,1) 
euro_garch1 <- garchFit(~ garch(1, 1), data = euro_res, trace = FALSE)
summary(euro_garch1)

gold_garch1 <- garchFit(~ garch(1, 1), data = gold_res, trace = FALSE)
summary(gold_garch1)

#Plot the fitted values against the observed values
euro_fit<- fitted(euro_model3)
gold_fit <- fitted(gold_model3)

par(mfrow=c(1,2))

plot(index(test_euro),test_euro, type = "l", col = "blue", lwd = 2, main = "Valeurs ajustées vs Observées - Euro/USD", ylab = "Retour", xlab = "Temps")
dates_fitted <- index(test_euro)[1:length(euro_fit)]
lines(dates_fitted, euro_fit, col = "red", lwd = 2)
legend("topright", legend = c("Valeurs Observées", "Valeurs Ajustées"), col = c("blue", "red"), lty = 1, lwd = 2)
grid()

plot(index(test_gold),test_gold, type = "l", col = "blue", lwd = 2, main = "Valeurs ajustées vs Observées - Gold/USD", ylab = "Retour", xlab = "Temps")
dates_fitted <- index(test_gold)[1:length(gold_fit)]
lines(dates_fitted, gold_fit, col = "red", lwd = 2)
legend("topright", legend = c("Valeurs Observées", "Valeurs Ajustées"), col = c("blue", "red"), lty = 1, lwd = 2)
grid()
#Check the residuals of the GARCH model

euro_garch_res <- residuals(euro_garch1, standardize = TRUE)
Box.test(euro_garch_res, lag = 20, type = "Ljung-Box")

gold_garch_res<- residuals(gold_garch1, standardize = TRUE)
Box.test(gold_garch_res, lag = 20, type = "Ljung-Box")

#Forcasting
#Forecast the next 20 observations

library(rugarch)

start <- as.Date("2024-09-01")
end <- as.Date("2024-09-20")

gold_20 <- getSymbols("GC=F", src = "yahoo", from = start, to = end, auto.assign = FALSE)
euro_20 <- getSymbols("EURUSD=X", src = "yahoo", from = start, to = end, auto.assign = FALSE)

gold_actual <- Cl(na.omit(gold_20))  # Closing prices for Gold/USD
euro_actual <- Cl(na.omit(euro_20)) 

# Spécification du modèle ARIMA-GARCH 
spec_euro <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(2, 3), include.mean = TRUE),
  distribution.model = "norm"
)

spec_gold <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),  
  mean.model = list(armaOrder = c(2, 2), include.mean = TRUE),
  distribution.model = "norm"  
)

# Ajustement du modèle ARIMA-GARCH 
euro_garch_fit <- ugarchfit(spec = spec_euro, data = euro_cl)
gold_garch_fit <- ugarchfit(spec = spec_gold, data = gold_cl)

# Prédiction des prochaines observations pour EUR/USD
euro_garch_forecast <- ugarchforecast(euro_garch_fit, n.ahead = 20)
gold_garch_forecast <- ugarchforecast(gold_garch_fit, n.ahead = 20)

# Convertir les prévisions en objets time series pour un meilleur affichage
euro_forecast_ts <- ts(as.numeric(fitted(euro_garch_forecast)), start = end(euro_cl), frequency = 1)
gold_forecast_ts <- ts(as.numeric(fitted(gold_garch_forecast)), start = end(gold_cl), frequency = 1)

#Intervalle de confiance 95%

euro_conf_int_lower <- euro_forecast_ts - 1.96 * sigma(euro_garch_forecast)
euro_conf_int_upper <- euro_forecast_ts + 1.96 * sigma(euro_garch_forecast)

gold_conf_int_lower <- gold_forecast_ts - 1.96 * sigma(gold_garch_forecast)
gold_conf_int_upper <- gold_forecast_ts + 1.96 * sigma(gold_garch_forecast)

# Plot Réel vs Prévisions
par(mfrow=c(1,1))

plot(index(euro_cl), euro_cl, main = "EURO/USD - Réel vs Prédit", col = "lightblue", xlab = "Temps", ylab = "Prix", type='l')
lines(index(euro_actual), euro_actual, col= "blue", lwd=2)
lines(euro_forecast_ts, col = "red", lwd=2)  # Ajouter les prévisionsyellow
legend("topleft", legend = c("Série Initiale", "Prédit","20 obs réelles"), col = c("lightblue", "red", "blue"), lty = 1)
grid()

par(mfrow=c(1,1))
plot(index(gold_cl), gold_cl, main = "GOLD/USD - Réel vs Prédit", col = "gold", xlab = "Temps", ylab = "Prix", type='l')
lines(index(gold_actual), gold_actual, col= "blue", lwd=2)
lines(gold_forecast_ts, col = "red", lwd=2)  # Ajouter les prévisions
legend("topleft", legend = c("Série Initiale", "Prédit","20 obs réelles"), col = c("gold", "red", "blue"), lty = 1)
grid()
#zoom

par(mfrow=c(1,1))
plot(index(euro_actual), euro_actual,ylim = c(min(euro_conf_int_lower), max(euro_conf_int_upper)), main = "EURO/USD - Réel vs Prédit", col = "blue", xlab = "Temps", ylab = "Prix", type='l')
lines(euro_forecast_ts, col = "red")  # Ajouter les prévisions
lines(euro_conf_int_lower, col = "green", lty = 2)
lines(euro_conf_int_upper, col = "green", lty = 2)
legend("topleft", legend = c("Réel", "Prédit","Bornes Inf/Sup"), col = c("blue", "red","green"), lty = 1)
grid()

par(mfrow=c(1,1))
plot(index(gold_actual), gold_actual, ylim = c(min(gold_conf_int_lower), max(gold_conf_int_upper)),main = "Gold/USD - Réel vs Prédit", col = "yellow", xlab = "Temps", ylab = "Prix", type='l')
lines(gold_forecast_ts, col = "red")  # Ajouter les prévisions
lines(gold_conf_int_lower, col = "green", lty = 2)
lines(gold_conf_int_upper, col = "green", lty = 2)
legend("topleft", legend = c("Réel", "Prédit","Bornes Inf/Sup"), col = c("yellow", "red","green"), lty = 1)
grid()

#Forecast Evaluation 1:
#Compare the forecasted values to the actual observations.

start <- as.Date("2024-09-01")
end <- as.Date("2024-09-20")

gold_20 <- getSymbols("GC=F", src = "yahoo", from = start, to = end, auto.assign = FALSE)
euro_20 <- getSymbols("EURUSD=X", src = "yahoo", from = start, to = end, auto.assign = FALSE)

gold_actual <- Cl(na.omit(gold_20))  # Closing prices for Gold/USD
euro_actual <- Cl(na.omit(euro_20)) 


# Mean Squared Error (MSE) and Mean Absolute Percentage Error (MAPE)
library(Metrics)

euro_actual_values <- as.numeric(euro_actual)
gold_actual_values <- as.numeric(gold_actual)

euro_forecast_values <- as.numeric(euro_forecast_ts)
gold_forecast_values <- as.numeric(gold_forecast_ts)

euro_mse <- mse(euro_actual_values, euro_forecast_values)
euro_mape <- mape(euro_actual_values, euro_forecast_values)

gold_mse <- mse(gold_actual_values, gold_forecast_values)
gold_mape <- mape(gold_actual_values, gold_forecast_values)

cat("EUR/USD Forecast Performance:\n")
cat("MSE:", euro_mse, "\n")
cat("MAPE:", euro_mape * 100, "%\n\n")  

cat("Gold/USD Forecast Performance:\n")
cat("MSE:", gold_mse, "\n")
cat("MAPE:", gold_mape * 100, "%\n\n") 

#scénario analysis


#alpha1 modifié 
spec_euro_scenario <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 2)), 
  mean.model = list(armaOrder = c(2, 3), include.mean = TRUE),
  #fixed.pars = list(alpha1 = 0.0367, beta1 = 0.98302),
  distribution.model = "norm"
)

spec_gold_scenario <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 2)), 
  mean.model = list(armaOrder = c(2, 2), include.mean = TRUE),
  #fixed.pars = list(alpha1 = 0.0367, beta1 = 0.94323),
  distribution.model = "norm"
)

# Ajustement du modèle ARIMA-GARCH avec nouveaux paramètres

euro_garch_fit_scenario <- ugarchfit(spec = spec_euro_scenario, data = euro_cl)
gold_garch_fit_scenario <- ugarchfit(spec = spec_euro_scenario, data = gold_cl)

euro_garch_fit_scenario <- ugarchfit(spec = spec_euro_scenario, data = euro_cl,
                                     #solver.control = list(trace = 0), 
                                     #fixed.pars = list(alpha1 = 0.148, beta1 = 0.98302)
                                     )
gold_garch_fit_scenario <- ugarchfit(spec = spec_euro_scenario, data = gold_cl,
                                     #solver.control = list(trace = 0), 
                                     #fixed.pars = list(alpha1 = 0.0367, beta1 = 0.94323)
                                     )

# Prédiction des prochaines observations (20 jours) dans les deux scénarios
euro_forecast_scenario <- ugarchforecast(euro_garch_fit_scenario, n.ahead = 20)
gold_forecast_scenario <- ugarchforecast(gold_garch_fit_scenario, n.ahead = 20)

# Extraire les prévisions et les convertir en séries temporelles
euro_forecast_scenario_ts <- ts(as.numeric(fitted(euro_forecast_scenario)), start = end(euro_cl), frequency = 1)
gold_forecast_scenario_ts <- ts(as.numeric(fitted(gold_forecast_scenario)), start = end(gold_cl), frequency = 1)

# Plot comparatif des scénarios pour EUR/USD
par(mfrow = c(2, 1))

plot(index(euro_actual), euro_actual, type = "l", col = "blue", lwd = 2, main = "EUR/USD - ORDRE = 2")
lines(ts(euro_forecast_scenario_ts, start = end(euro_cl)), col = "red", lwd = 2)
legend("topleft", legend = c("Réel", "Prévision"), col = c("blue", "red"), lty = 1, lwd = 2)
grid()

plot(index(euro_actual), euro_actual, type = "l", col = "blue", lwd = 2, main = "EUR/USD - ORDRE = 1")
lines(ts(euro_forecast_ts, start = end(euro_cl)), col = "red", lwd = 2)
legend("topleft", legend = c("Réel", "Prévision ordre Modifié"), col = c("blue", "red"), lty = 1, lwd = 2)
grid()

# Plot comparatif des scénarios pour Gold/USD
par(mfrow = c(2, 1))

plot(index(gold_actual), gold_actual, type = "l", col = "gold", lwd = 2, main = "Gold/USD - Garch(1,2)")
lines(ts(gold_forecast_scenario_ts, start = end(gold_cl)), col = "red", lwd = 2)
legend("topleft", legend = c("Réel", "Prévision"), col = c("gold", "red"), lty = 1, lwd = 2)
grid()

plot(index(gold_actual), gold_actual, type = "l", col = "gold", lwd = 2, main = "Gold/USD - Garch(1,1)")
lines(ts(gold_forecast_ts, start = end(gold_cl)), col = "red", lwd = 2)
legend("topleft", legend = c("Réel", "Prévision ordre modifié"), col = c("gold", "red"), lty = 1, lwd = 2)
grid()

#Forecasting Evaluation 2
#BOUCLE FORECAST with ARIMA only

# Étendre la série de données avec les prévisions pas à pas
euro_cl_extended <- euro_cl
gold_cl_extended <- gold_cl

# Déterminer les dates futures (au-delà de la dernière date de votre jeu de données)
start_date1 <- index(euro_cl)[length(euro_cl)] + 1  # Date juste après la dernière observation
future_dates1 <- seq(start_date1, by = "days", length.out = 20)  # 20 prévisions, par exemple, quotidiennes

start_date2 <- index(gold_cl)[length(gold_cl)] + 1  # Date juste après la dernière observation
future_dates2 <- seq(start_date2, by = "days", length.out = 20)  # 20 prévisions, par exemple, quotidiennes

# Boucle pour générer les prévisions en ajoutant chaque valeur au modèle
for (i in 1:20) {
  # Générer la prévision pour l'étape i
  euro_forecast <- forecast(euro_model2, h = 1, level = 95)
  euro_forecast_values <- euro_forecast$mean
  
  # Ajouter la prévision à la série étendue avec la date correspondante
  euro_cl_extended <- rbind(euro_cl_extended, xts(euro_forecast_values, order.by = future_dates1[i]))
  
  # Mettre à jour le modèle pour intégrer la nouvelle prévision (réajustement facultatif)
  euro_model2 <- Arima(euro_cl_extended, order = c(2, 1, 3), include.mean = TRUE)
}

# Boucle GOLD
for (i in 1:20) {
  # Générer la prévision pour l'étape i
  gold_forecast <- forecast(gold_model2, h = 1, level = 95)
  gold_forecast_values <- gold_forecast$mean
  
  # Ajouter la prévision à la série étendue avec la date correspondante
  gold_cl_extended <- rbind(gold_cl_extended, xts(gold_forecast_values, order.by = future_dates2[i]))
  
  # Mettre à jour le modèle pour intégrer la nouvelle prévision (réajustement facultatif)
  gold_model2 <- Arima(gold_cl_extended, order = c(2, 1, 2), include.mean = TRUE)
}

# Afficher les prévisions et la série étendue
print(euro_cl_extended)
print(gold_cl_extended)

euro_forecast_20 <- tail(euro_cl_extended, 20)
gold_forecast_20 <- tail(gold_cl_extended, 20)

# Visualisation des prévisions

par(mfrow=c(1,1))
plot(index(euro_cl), euro_cl, type = "l", col = "lightblue", lwd = 2, main = "Prévisions ARIMA - Euro/USD", ylab = "Taux de change", xlab = "Temps")
lines(index(euro_actual), euro_actual, col= "blue", lwd=2)
lines(index(euro_forecast_20), euro_forecast_20, col = "red", lwd = 2)  # Afficher les prévisions sur la série étendue
legend("topleft", legend = c("Série initiale","20 obs réelles", "Prévisions ARIMA"), col = c("lightblue","blue", "red"), lty = 1, lwd = 2)
grid()

par(mfrow=c(1,1))
plot(index(gold_cl), gold_cl, type = "l", col = "gold", lwd = 2, main = "Prévisions ARIMA - Gold/USD", ylab = "Taux de change", xlab = "Temps")
lines(index(gold_actual), gold_actual, col= "blue", lwd=2)
lines(index(gold_forecast_20), gold_forecast_20, col = "red", lwd = 2)  # Afficher les prévisions sur la série étendue
legend("topleft", legend = c("Série initiale","20 obs réelles", "Prévisions ARIMA"), col = c("gold","blue", "red"), lty = 1, lwd = 2)
grid()

#zoom
par(mfrow=c(1,1))
plot(index(euro_actual), euro_actual, type = "l", col = "blue", lwd = 2, main = "Prévisions ARIMA - Euro/USD", ylab = "Taux de change", xlab = "Temps")
lines(index(euro_forecast_20), euro_forecast_20, col = "red", lwd = 2)  # Afficher les prévisions sur la série étendue
legend("topleft", legend = c("Observations réelles", "Prévisions ARIMA"), col = c("blue", "red"), lty = 1, lwd = 2)
grid()

par(mfrow=c(1,1))
plot(index(gold_actual), gold_actual, type = "l", col = "gold", lwd = 2, main = "Prévisions ARIMA - Gold/USD", ylab = "Taux de change", xlab = "Temps")
lines(index(gold_forecast_20), gold_forecast_20, col = "red", lwd = 2)  # Afficher les prévisions sur la série étendue
legend("topleft", legend = c("Observations réelles", "Prévisions ARIMA"), col = c("gold", "red"), lty = 1, lwd = 2)
grid()


#mse forecast ARIMA onlY

euro_mse2 <- mse(euro_actual, euro_forecast_20)
euro_mape2 <- mape(euro_actual, euro_forecast_20)

gold_mse2 <- mse(gold_actual, gold_forecast_20)
gold_mape2 <- mape(gold_actual, gold_forecast_20)

cat("EUR/USD Forecast Performance:\n")
cat("MSE:", euro_mse2, "\n")
cat("MAPE:", euro_mape2 * 100, "%\n\n")  

cat("Gold/USD Forecast Performance:\n")
cat("MSE:", gold_mse2, "\n")
cat("MAPE:", gold_mape2 * 100, "%\n\n") 

# Visualisation des 2 PREV

start <- as.Date("2022-01-01")
end <- as.Date("2024-09-20")

gold_sum <- getSymbols("GC=F", src = "yahoo", from = start, to = end, auto.assign = FALSE)
euro_sum <- getSymbols("EURUSD=X", src = "yahoo", from = start, to = end, auto.assign = FALSE)

gold_sum <- Cl(na.omit(gold_sum))  # Closing prices for Gold/USD
euro_sum <- Cl(na.omit(euro_sum))

par(mfrow=c(1,1))
plot(index(euro_sum), euro_sum, type = "l", col = "blue", lwd = 2, main = "Comparaison des Prévisions - Euro/USD", ylab = "Taux de change", xlab = "Temps")
lines(index(euro_forecast_20), euro_forecast_20, col = "red", lwd = 2)  # Afficher les prévisions sur la série étendue
lines(euro_forecast_ts, col = "green", lwd=2)
legend("topleft", legend = c("Observations réelles", "Prévisions ARIMA","Prévision ARIMA-GARCH"), col = c("blue", "red","green"), lty = 1, lwd = 2)
grid()

par(mfrow=c(1,1))
plot(index(gold_sum), gold_sum, type = "l", col = "gold", lwd = 2, main = "Comparaison des Prévisions - Gold/USD", ylab = "Taux de change", xlab = "Temps")
lines(index(gold_forecast_20), gold_forecast_20, col = "red", lwd = 2)  # Afficher les prévisions sur la série étendue
lines(gold_forecast_ts, col = "green", lwd=2)
legend("topleft", legend = c("Observations réelles", "Prévisions ARIMA", "Prévision ARIMA-GARCH" ), col = c("gold", "red","green"), lty = 1, lwd = 2)
grid()
#zoom

par(mfrow=c(1,1))
plot(index(euro_actual), euro_actual, type = "l", col = "blue", lwd = 2, main = "Comparaison des Prévisions - Euro/USD", ylab = "Taux de change", xlab = "Temps")
lines(index(euro_forecast_20), euro_forecast_20, col = "red", lwd = 2)  # Afficher les prévisions sur la série étendue
lines(euro_forecast_ts, col = "green", lwd=2)
legend("topleft", legend = c("Observations réelles", "Prévisions ARIMA","Prévision ARIMA-GARCH"), col = c("blue", "red","green"), lty = 1, lwd = 2)
grid()

par(mfrow=c(1,1))
plot(index(gold_actual), gold_actual, type = "l", col = "gold", lwd = 2, main = "Comparaison des Prévisions - Gold/USD", ylab = "Taux de change", xlab = "Temps")
lines(index(gold_forecast_20), gold_forecast_20, col = "red", lwd = 2)  # Afficher les prévisions sur la série étendue
lines(gold_forecast_ts, col = "green", lwd=2)
legend("topleft", legend = c("Observations réelles", "Prévisions ARIMA", "Prévision ARIMA-GARCH" ), col = c("gold", "red","green"), lty = 1, lwd = 2)
grid()


