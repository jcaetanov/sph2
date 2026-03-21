# SPH2 – Classical Tube Problem

## Abstract

Este repositório apresenta uma implementação em C++ do método Smoothed Particle Hydrodynamics (SPH) aplicada a um problema clássico de dinâmica dos fluidos: o escoamento em tubo com condições iniciais descontínuas.

O objetivo do projeto é validar numericamente o método SPH por meio da simulação de fenômenos como ondas de choque e rarefação, amplamente utilizados como benchmarks em hidrodinâmica computacional.

---

## 1. Introdução

Problemas unidimensionais em geometria de tubo são amplamente utilizados na validação de métodos numéricos para dinâmica dos fluidos. Entre eles, destacam-se os problemas de tubo de choque, que envolvem descontinuidades iniciais em propriedades como pressão e densidade.

Este projeto implementa uma solução baseada em SPH, um método Lagrangiano e sem malha, no qual o fluido é representado por partículas que interagem localmente.

---

## 2. Fundamentação Teórica

O método Smoothed Particle Hydrodynamics (SPH) é uma técnica numérica mesh-free que discretiza o fluido em partículas.

Cada propriedade física é aproximada por:

- Interpolação via funções kernel
- Soma das contribuições das partículas vizinhas

As equações resolvidas incluem:

- Conservação de massa
- Conservação de momento
- Equação de estado

Este tipo de abordagem é particularmente eficaz em problemas com grandes deformações e descontinuidades.

---

## 3. Problema do Tubo Clássico

O problema consiste em:

- Um domínio unidimensional representando um tubo
- Duas regiões com diferentes condições iniciais
- Evolução temporal que gera:
  - ondas de choque
  - ondas de rarefação
  - descontinuidades de contato

Este problema é amplamente utilizado para:

- Verificação de estabilidade numérica
- Avaliação da precisão do método
- Comparação com soluções analíticas

---

## 4. Estrutura do Código

### 4.1 Inicialização

O código define:

- Distribuição inicial das partículas
- Parâmetros físicos (densidade, pressão, energia)
- Condições iniciais distintas ao longo do domínio

---

### 4.2 Representação das Partículas

Cada partícula armazena:

- Posição
- Velocidade
- Densidade
- Pressão
- Energia interna

---

### 4.3 Kernel de Suavização

O método utiliza uma função kernel para interpolação das propriedades físicas.

Funções típicas incluem:

- Cubic spline
- Gaussian kernel

---

### 4.4 Cálculo de Densidade

A densidade é obtida por:

- Soma ponderada das partículas vizinhas
- Dependência do raio de suavização (smoothing length)

---

### 4.5 Cálculo de Forças

As forças consideradas incluem:

- Gradiente de pressão
- Termos viscosos (quando aplicável)

---

### 4.6 Integração Temporal

A evolução do sistema é realizada por métodos explícitos, como:

- Euler
- Leapfrog

---

### 4.7 Condições de Contorno

O domínio do tubo pode incluir:

- Fronteiras fixas
- Condições refletivas

---

### 4.8 Saída de Dados

O código gera dados para análise, como:

- Perfis de densidade
- Velocidade
- Pressão ao longo do tempo

---

## 5. Metodologia

O desenvolvimento seguiu as etapas:

1. Definição do problema físico
2. Discretização do domínio em partículas
3. Implementação das equações do SPH
4. Simulação temporal
5. Comparação qualitativa com soluções conhecidas

---

## 6. Resultados e Discussão

A simulação permite observar:

- Propagação de ondas de choque
- Formação de regiões de rarefação
- Sensibilidade a parâmetros numéricos

Resultados típicos incluem:

- Perfis não lineares de densidade e pressão
- Estruturas características de soluções hiperbólicas

---

## 7. Limitações

- Sensibilidade à escolha do kernel
- Dependência do passo de tempo
- Custo computacional elevado para alta resolução

---

## 8. Conclusão

A implementação demonstra a capacidade do método SPH em resolver problemas clássicos de dinâmica dos fluidos, validando sua aplicação em cenários com descontinuidades e grandes gradientes.

Este tipo de benchmark é essencial para o desenvolvimento de simuladores mais complexos e para aplicações em pesquisa científica.

---

## 9. Execução

### Compilação
