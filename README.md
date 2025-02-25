"North Star" rendering reference project
----------------------------------------

The purpose of this project is to provide the basic skeleton of the rendering
capabilities LANL is looking for so that we can identify parts of `ggmolvis`
that need improvement. Generally, these pain points will be indicated with `# NOTE:`
in the rendering source code here--for example, when we need to manually use `bpy`
to accomplish something we cannot currently with `ggmolvis`. Then we'll typically
aim to open an upstream issue or PR to close the capability gap and make the
ggmolvis rendering scripts easier to write/understand long-term.

We use a YAML input deck format here to specify various rendering options and
the paths to the trajectory files so that in principle we can use this repo
to point to rendering challenges for both publicly-available trajectories
and in-house trajectories without changing the Python source (just the input
deck).

Sample incantation on command line:

``
