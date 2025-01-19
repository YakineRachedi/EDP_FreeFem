# Résolution d'un problème de transfert de chaleur :

Soit $\Omega$ un domaine borné suffisamment régulier, on considère le problème :

$$
\begin{cases}
    -k \Delta \phi + \vec{u} \cdot \nabla \phi = 0 & \text{dans } \Omega, \\
    \phi = 0 & \text{sur } \Gamma_c, \\
    \phi = 1 & \text{sur } \Gamma_H, \\
    \nabla \phi \cdot \vec{n} = 0 & \text{sur } \Gamma_w.
\end{cases}
$$

Le domaine $\Omega$ est illustré comme suit :
![Illustration du domaine $\Omega$](GeometrieOmega.PNG)

<figcaption style="text-align: center; font-weight: bold;">
Figure 1 – Illustration du domaine $\Omega$.
</figcaption>


Illustration du domaine ($\Omega$)

# Formulation Variationnelle

Soit $\Omega$ un domaine borné suffisamment régulier, on cherche $\phi \in H^1(\Omega)$ tel que :

1. $\phi = 0$ sur $\Gamma_c$ (condition de Dirichlet).
2. Pour tout $v \in H^1_0(\Omega)$, on a l'égalité suivante :

$$
\int_\Omega k \nabla \phi \cdot \nabla v \, dx + \int_\Omega (\vec{u} \cdot \nabla \phi) v \, dx = 0
$$

où :

- $H^1(\Omega)$ est l'espace de Sobolev.
- $H^1_0(\Omega)$ est le sous-espace de $H^1(\Omega)$ des fonctions nulles sur $\Gamma_c$.

Les conditions aux limites sont intégrées dans les espaces fonctionnels ou directement dans l'intégrale.

- $\Gamma_H$ impose $\phi = 1$ (on peut ajuster ce terme en tenant compte des conditions de Dirichlet avec des fonctions de test modifiées).
- $\Gamma_w$ impose $\nabla \phi \cdot \vec{n} = 0$ et est prise en compte naturellement dans le cadre variationnel, car elle correspond à une condition de Neumann homogène.

D'où :

$$
a(\phi, v) = \int_\Omega k \nabla \phi \cdot \nabla v \, dx + \int_\Omega (\vec{u} \cdot \nabla \phi) v \, dx
$$

et

$$
l(v) = 0.
$$

# Problème des Éléments Finis

Pour $ \Gamma_c = \Gamma_f $, $ \Gamma_H = \Gamma_2 $, et $ \Gamma_w = \Gamma_1 \cup \Gamma_3 $ :

- Cas 1 : $k = 0.01$ et $\vec{u} = \vec{0}$.
- Cas 2 : $k = 0.01$ et $\vec{u} = (0.1, 0)$.

Nous cherchons à résoudre le problème défini précédemment dans **FreeFem++** en utilisant des éléments finis **P1** sur l'un des maillages, par exemple, le maillage 4 avec $ h = 0.1 $.

# Étude de convergence :

Pour le **cas 1**, nous allons étudier la convergence. On note $\phi_{\text{ref}}$ la solution de référence du problème précédent sur le maillage 1, et $\phi_2$, $\phi_3$, $\phi_4$ les solutions obtenues respectivement sur les maillages 2, 3 et 4.

Les erreurs absolues et relatives, pour $ i = 2, 3, 4 $, sont définies par :

$$
E_i^A = \phi_{\text{ref}} - \Pi_1^i \phi_i \quad \text{et} \quad E_i^R = \frac{\| E_i^A \|}{\| \phi_{\text{ref}} \|_{L^2(\Omega)}}
$$

où $\Pi_1^i$ est l'opérateur d'interpolation de l'espace des éléments finis associé au maillage $i$ vers l'espace des éléments finis associé au maillage 1.
