
This is a log for tests of the single scattering equation

I want to create a tests that proves that something might be wrong there or in the multiple scattering 


Solver:
_______

Call graph:

- single_scatterer_nonrel
- semrel (There is even some potfit going on)
- rwave (This called burlisch-stoer steper)


Here the problem starts

The main goal would be to replace semrel. This routine seems to be flawed.





Tests carried out:
------------------

1. remove XC: [x]
2. remove vmt: Reduces the problem

My suggestion is that the solver is wrong. This is amplifying all issue

3. Code is converging if one uses broyden of potential
4. Core states are absolutly accurate. They are also normalized so.
5. Fix core states: [x] also doesn't help
