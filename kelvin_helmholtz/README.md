# Kelvin-Helmholtz instability simulation.

## pmocz的代码

[github仓库地址](https://github.com/pmocz/KelvinHelmholtzInstability/blob/master/KelvinHelmholtzInstability.py)

73-184行, 满屏的公式, Pmocz大哥真是有毅力...

用python2.7写的. 居然能跑通, 太恐怖了... 简直写得和汇编一样, 让人心生惧意.

## Kelvin-Helmholtz 不稳定性的原理

## comsol

初始速度值:

$$
x=
tanh(30*(y[1/m]-0.25))*(y<=0.5)+tanh(30*(0.75-y[1/m]))*(y>0.5)
$$

$$
y=
0.05*sin(2*pi*x[1/m])
$$

