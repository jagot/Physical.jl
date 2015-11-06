using Physical, Quantities, Base.Test
using Physical.NaturalUnits

au = get_natural_unity(:Hartree)
auRy = get_natural_unity(:Rydberg)
planck = get_natural_unity(:Planck)
# natural = get_natural_unity(:Natural)
# geom = get_natural_unity(:Geometrized)

function test_conv(v,s=au)
    v_s = s(v)
    println(v_s.value ≈ 1.0)
end

function test_exact_unity(v,s=au)
    v_s = s(v)
    @test v_s.value == 1
end

function test_approx_unity(v, ε = 1e-7, s=au)
    v_s = s(v)
    @test_approx_eq_eps v_s.value 1.0 ε
end

test_exact_unity(Hbar_planck)
test_exact_unity(e_electron)
test_exact_unity(m_electron)
test_exact_unity(K_coulomb)
test_exact_unity(K_boltzmann)
test_approx_unity(A_bohr, 1e-9)
test_approx_unity(2.18769126e6Meter/Second, 1e-8)
test_approx_unity(1.99285166e-24Kilogram*Meter/Second)
test_approx_unity(2.41888430e-17Second)
test_approx_unity(4.13413738e16Hertz)
test_approx_unity(4.3597442e-18Joule)
test_approx_unity(27.211ElectronVolt, 1e-4)
test_approx_unity(E_hartree)
test_approx_unity(5.14220651e11Volt/Meter)

lP = √(Hbar_planck*G_gravitation/C_light^3)
test_exact_unity(lP, planck)
