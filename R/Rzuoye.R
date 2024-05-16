#' ODE solution
#'
#' This function solves an ODE system based on the given parameter list.
#' @param user created parameters
#' @return solution of an ODE function system
#' @export
library(deSolve)
library(ggplot2)
# 扩展的ODE模型函数
model_extended <- function(t, y, param) {
  S <- y[1]
  E <- y[2]
  I1 <- y[3]
  I2 <- y[4]
  R <- y[5]
  N <- param$N
  beta1 <- param$beta1
  beta2 <- param$beta2
  mu <- param$mu
  gamma1 <- param$gamma1
  gamma2 <- param$gamma2
  lamda <- param$lamda
  dSt <- mu * (N - S) - beta1 * S * I1 / N - beta2 * S * I2 / N
  dEt <- beta1 * S * I1 / N + beta2 * S * I2 / N - mu * E - lamda * E
  dI1t <- lamda * E - mu * I1 - gamma1 * I1
  dI2t <- gamma1 * I1 - mu * I2 - gamma2 * I2
  dRt <- gamma1 * I1 + gamma2 * I2 - mu * R
  outcome <- c(dSt, dEt, dI1t, dI2t, dRt)
  return(list(outcome))
}
# 时间序列
times <- seq(0, 156, by = 1/7)

# 参数设置
param <- list(
  N = 1,
  beta1 = 4,
  beta2 = 3, # 假设二次感染的传播率略低
  mu = 0.000,
  gamma1 = 0.1,
  gamma2 = 0.08, # 假设二次感染的恢复率略低
  lamda = 0.03
)

# 初始状态
init <- c(S = 0.9999, E = 0.00008, I1 = 0.00001, I2 = 0.00001, R = 0)

# 解决ODE系统
result <- ode(y = init, times = times, func = model_extended, parms = param)

# 转换为数据框
result <- as.data.frame(result)

# 绘制结果
seirplot_extended <- ggplot(data = result) +
  geom_line(aes(x = time, y = S, col = 'S'), lwd = 2) +
  geom_line(aes(x = time, y = E, col = 'E'), lwd = 2) +
  geom_line(aes(x = time, y = I1, col = 'I1'), lwd = 2) +
  geom_line(aes(x = time, y = I2, col = 'I2'), lwd = 2) +
  geom_line(aes(x = time, y = R, col = 'R'), lwd = 2) +
  labs(x = 'Time', y = 'Ratio') +
  scale_color_manual(name = 'SEIR Extended',
                     values = c("S" = "orange", "E" = "purple", "I1" = "blue", "I2" = "cyan", "R" = "green"))

# 显示图形
print(seirplot_extended)

# 保存图形为PDF
ggsave(seirplot_extended, file = 'lubing_extended.pdf', width = 7, height = 6)

# 保存图形为SVG
ggsave(seirplot_extended, file = 'lubing_extended.svg', width = 7, height = 6)




