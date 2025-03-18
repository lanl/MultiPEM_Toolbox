Multi-Phenomenology Explosion Monitoring (MultiPEM) Toolbox

The MultiPEM approach involves the statistical modeling of observed
signatures from multiple phenomenologies at distributed sensor networks
to infer characteristics of interest from new explosive events. Recent
publications have considered processed waveforms from seismic (e.g.
P-wave displacement, peak particle velocity), acoustic (e.g. pressure
impulse, duration), and optical (e.g. irradiance first minimum, second
maximum) phenomenologies for this purpose. Expected signatures are
modeled empirically or from first principles utilizing two basic
components: source physics which generates initial waveforms from
explosive events, and path physics which attenuates these waveforms
along their propagation route from source to sensor. Systematic
departures of observed signatures from these two model components can
be modelled dynamically.

Historically, device parameters were characterized by single
phenomenology analyses. When signatures from multiple phenomenologies
are available, MultiPEM assessments will naturally compute a compromise
among single phenomenology results while often producing rigorously
quantified uncertainty reduction as the desired benefit.

Consult the technical report provided with this repository for a
description of the statistical methodology and example application
(see ./Runfiles/IYDT) implemented in this repository:

./multipem-la-ur-23-21950.pdf

Consult the user manual provided with this repository for descriptions
of the input decks provided with the example application of this
technical report (see ./Runfiles/IYDT), and more generally instructions
on how to run single-phenomenology and MultiPEM analyses:

./multipem-um-la-ur-23-30117.pdf

Consult the presentation prepared for the AFTAC Geophysical Sciences
Review Panel for extensions of "nested" (within source) path effects
modeling to "crossed" (across source) path effects modeling also
implemented in this repository (see ./Runfiles/IYDT-gsrp):

./MultiPEM-GSRP-030425.pdf

Directories in this repository:

Applications/
	Contains application-specific code, data, and
	verification tests

Code/
	Contains functions utilized across applications,
	such as the log-likelihood, log-prior, and log-
	posterior distributions and their associated
	gradients

Runfiles/
	Set up application directories to run
	MultiPEM analyses using R

Runfiles-Docker/
	Set up application directories to run
	MultiPEM analyses using R from Docker

Test/
	Contains code for running general verification
	tests using R

Test-Docker/
	Contains code for running general verification
	tests using R from Docker

These directories, and associated sub-directories, contain additional
README files with information and instructions for running the code.

This toolbox is open source under the BSD-3 License.
Its LANL-internal identifier is O4673.
