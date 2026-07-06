---
name: New Project Proposal
about: This form is for new software projects that like to join HPSF.
title: "Proposal for LAMMPS to join the High Performance Software Foundation,"
labels: [Project Proposal]
author: The LAMMPS Development Team, developers@lammps.org
assignees: ''
urlcolor: blue
---


### 1. Name of Project

LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

### 2. Project Description

LAMMPS is a particle simulation code with a focus on large-scale
classical molecular dynamics for materials.  It is widely used in
materials science, materials chemistry, condensed matter physics, and
related fields.  LAMMPS is designed to work well on massively parallel
computers, while also achieving good single-node performance on CPU and
GPU systems.  It has a large, well-established community of users
throughout the world.  LAMMPS development is coordinated through a
formal process of pull requests, code review, testing, and releases
managed via GitHub.

LAMMPS is deployed and used in production for academic research at all
major DOE computing facilities — including the leadership computing
facilities at Oak Ridge and Argonne, as well as Livermore, Los Alamos,
Sandia, and NERSC — and at NSF-sponsored supercomputing centers like
TACC and PSC.  Adequate LAMMPS performance is frequently among the
acceptance criteria for new hardware deployments at such centers.
Hardware vendors including NVIDIA, AMD, and Intel regularly dedicate
software engineering staff time to contribute performance enhancements
and improved support for their hardware, and HPC center staff have
contributed improvements that lift parallel scaling limitations.  LAMMPS
is also used in industrial research.

The diversity of participation in LAMMPS *development* is documented by
its list of contributing authors.  The project's authors page,
[authors.html](https://www.lammps.org/authors.html),
records more than two decades of
feature contributions in a continuous timeline: over 200 recognized
contributors of significant features, drawn from national laboratories,
universities, HPC centers, and companies around the world.  More
recently, the provenance of significant contributions is publicly
documented in the release notes at
[lammps/releases](https://github.com/lammps/lammps/releases).
The reach of LAMMPS extends
well beyond its developer community, however: over the years its
audience has shifted from users with software development experience
toward non-technical users who rely on HPC center installations or
pre-compiled binary packages.  The size of this wider user base is
difficult to quantify directly; the citation statistics below are its
most visible indicator.

The original 1995 LAMMPS article was cited about 5,100 times in 2021
alone; since the publication of the second canonical article
[(Comp Phys Comm, 2022)](https://doi.org/10.1016/j.cpc.2021.108171), the two articles combined draw
roughly 8,000 citations per year.  There were approximately 405,000
downloads of LAMMPS source tarballs from the project website between
September 2004 and mid-2021.  Since distribution moved to GitHub — where
source downloads are no longer individually countable — the pre-compiled
binary and LAMMPS-GUI packages alone (for Windows, macOS, and Linux)
have been downloaded about 26,000 times across the last 33 releases.
LAMMPS is included in the Spack package manager (where it is maintained
by a LAMMPS core developer who also serves on the Spack management
team), in conda-forge, and in major Linux distributions.

### 3. Statement on Alignment with High Performance Software Foundation's Mission

Since its creation three decades ago, LAMMPS development has been
concentrated at Sandia National Laboratories' Center for Computing
Research and for nearly two decades at Temple University's Institute for
Computational Molecular Science.  LAMMPS now has core developers and
important stakeholders at multiple additional institutions in the United
States.  In order to ensure the long-term health of LAMMPS, we wish to
establish a governance model that is not dependent on any one person or
institution.

Given LAMMPS's history, sustained growth, broad base of
multi-institutional contributors and users, and its role as widely
adopted, production-grade infrastructure across academia, national
laboratories, and industry, we are requesting consideration for
admission at the **Core Stage**.

### 4. Project Website

[lammps.org](https://www.lammps.org)

### 5. Open Source License

LAMMPS is distributed under the GNU General Public License version 2
([SPDX: GPL-2.0-only](https://spdx.org/licenses/GPL-2.0-only.html)); see
[LICENSE](https://github.com/lammps/lammps/blob/develop/LICENSE).  A
variant under the GNU Lesser General Public License version 2.1, with
the small number of non-conforming parts removed, is available on
request.

### 6. Code of Conduct

[Code of Conduct](https://github.com/lammps/lammps/blob/develop/.github/CODE_OF_CONDUCT.md)
(adapted from the Contributor Covenant v1.4; reports are handled
confidentially via developers@lammps.org)

### 7. Governance Practices

LAMMPS development practices — the pull request, review, and testing
workflow, and the roles of core developers — are documented in
[`.github/CONTRIBUTING.md`](https://github.com/lammps/lammps/blob/develop/.github/CONTRIBUTING.md)
in the LAMMPS repository and also in the [Development
Processes](https://docs.lammps.org/Developer_org.html) sections of the
LAMMPS manual.  As part of joining HPSF, the project is formalizing its
long-standing informal governance practices in a standalone governance
document (defining the Contributor, Core Developer, and Lead Developer
roles, committer onboarding and emeritus processes, and consensus-based
decision-making); it will be published as `.github/GOVERNANCE.md` in the LAMMPS
repository.  A complete draft accompanies this proposal.

### 8. Two Sponsors from the High Performance Software Foundation's Technical Advisory Committee

- Curt Ober
- Andrew Myers

### 9. What is the project's solution for source control?

Git, hosted on GitHub: [github.com/lammps/lammps](https://github.com/lammps/lammps).  The `develop`
branch is protected: all changes must be submitted as pull requests,
must pass automated testing (GitHub Actions plus the project's own
continuous integration server at [ci.lammps.org](https://ci.lammps.org)), and must
receive review and approval by a core developer before merging.  Direct
pushes are not possible.  The release branches (`release` from `develop`
and `stable` from `maintenance`) advance only by fast-forward merges as
described in Section 12.

### 10. What is the project's solution for issue tracking?

GitHub Issues on the main repository:
[github.com/lammps/lammps/issues](https://github.com/lammps/lammps/issues), used for bug reports and
feature requests, with discussion of proposed changes taking place on
the corresponding pull requests.

### 11. Please list all external dependencies and their license

The LAMMPS core has no hard external dependencies beyond a C++17
compiler; an MPI library (e.g., Open MPI, MPICH, Intel MPI) is recommended
for parallel execution but optional.  All further dependencies are tied to
individual optional packages (94 in the current release); the
authoritative per-package list of libraries, versions, and build instructions
is maintained in the manual at [Build extras](https://docs.lammps.org/Build_extras.html).

Representative examples:

| Library | Used by | License |
|---|---|---|
| FFTW3 | KSPACE (long-range electrostatics) | GPL-2.0-or-later |
| Eigen3 | several fitting/ML packages | MPL-2.0 |
| Kokkos (bundled) | KOKKOS performance-portable backend | Apache-2.0 with LLVM exception |
| Voro++ | VORONOI package | BSD-3-Clause |
| libpng / zlib / libjpeg | `dump image`, compressed I/O | libpng / zlib / IJG licenses |
| COLVARS (bundled) | COLVARS package (enhanced sampling) | LGPL-3.0 |
| PLUMED | PLUMED package (enhanced sampling) | LGPL-3.0 |
| KIM API | KIM package (interatomic potentials) | LGPL-2.1 |

When the CMake build system downloads an external library automatically,
the download is validated against a pinned SHA-256 checksum or Git tag.
For these downloads the project additionally maintains vetted copies on
its own download servers as a fallback, should the upstream source
become temporarily or permanently unavailable.

### 12. Please describe your release methodology and mechanics

LAMMPS follows a continuous release model built on four branches:

- All accepted changes are merged into the protected `develop` branch
  via reviewed, CI-passing pull requests.  The project aims to keep
  `develop` fully functional at all times; additional tests run after
  merging, and all tests must pass before any release is made.
- Every 4–8 weeks, `develop` is merged into the `release` branch with a
  fast-forward merge and tagged (`patch_<date>`) as a **feature release**.
- Approximately once per year, integrating new features is deferred in
  favor of extended testing, fixing bugs, and reviewing documentation
  followed by a **stable release**: the `stable` branch is then reset
  to the corresponding state of `release` and receives an additional
  `stable_<date>` tag.
- Between stable releases, bug fixes are backported from `develop` to the
  `maintenance` branch and periodically fast-forward-merged into `stable`,
  published as **update releases** (`stable_<date>_update<N>`).  The
  `maintenance` branch is reset from `release` or `stable` at each new
  stable release.

The full release procedure is documented in
[`.github/release_steps.md`](https://github.com/lammps/lammps/blob/develop/.github/release_steps.md)
in the LAMMPS repository.
Since the 10 September 2025 release,
GitHub releases are configured as immutable: neither the tag nor the
published assets can be altered after release, and GitHub generates a
release attestation that can be used to verify asset integrity.  SHA-256
checksums are additionally published for all files downloadable from
[lammps.org](https://www.lammps.org).

### 13. Please describe Software Quality efforts (CI, security, auditing)

- **Continuous integration:** every pull request must pass an extensive
  matrix of automated checks before merging, using GitHub Actions
  workflows (compilation and unit tests on Linux, macOS, and
  MSVC/Windows; compatibility check with 64-bit atom indexing, single
  precision FFTs, the legacy build system, the C++23 standard; and
  checks for documented programming-style conformance).  Additional,
  more time consuming tests are run after a merge and required to be
  cleared for a release (KOKKOS builds with single/mixed/double
  floating-point precision; MSVC/Windows compilation;) Compilations use
  high warning levels (`-Wall -Wextra`) with GCC and Clang.
- **Testing:** unit tests are implemented with GoogleTest (currently 615–769
  tests per configuration across about a dozen platform/configuration
  combinations).  A physics-level regression test suite based on the bundled
  example inputs runs via dedicated quick- and full-regression workflows;
  this suite is a recent addition and expanding its coverage is ongoing work.
- **Testing dashboard:** current results for all of the above —
  including unit test status, regression test progress, and code
  coverage — are published on the LAMMPS Testing Dashboard at
  [lammps-test-results](https://lammps.github.io/lammps-test-results/)
- **Static analysis:** CodeQL (with the security-and-quality query
  suite) runs on the C++ and Python code in CI; Coverity Scan, the Clang
  Static Analyzer, and clang-tidy run regularly on dedicated
  infrastructure, with results tracked on the testing dashboard.  Since
  systematic static analysis was adopted in 2020, roughly 5,000 reported
  defects have been reviewed and fixed, and the outstanding-defect count
  has been driven down from over 2,300 to under 500.
- **Code review:** every change requires approval by a core developer
  other than its author, and merging is performed by a third party (a
  core developer other than the approver), so at least two core
  developers see every contribution.
- **Security:** the project security policy is published at
  [SECURITY.md](https://github.com/lammps/lammps/blob/develop/SECURITY.md).
  It documents the threat model appropriate for user-level scientific
  software, the reporting path for suspected vulnerabilities, and the
  supply-chain integrity measures described in Sections 11 and 12 (checksum
  validation of dependency downloads, SHA-256 checksums for all published
  archives, immutable GitHub releases with attestation).  To guard against
  source obfuscation attacks using Unicode homoglyphs or bidirectional
  reordering, only all-ASCII source code is accepted; this is enforced by
  automated checks on every submission.

### 14. Please list the project's leadership team

- Steve Plimpton, Sandia National Laboratories (retired)
- Axel Kohlmeyer, Temple University
- Aidan Thompson, Sandia National Laboratories

### 15. Please list the project members with access to commit to the mainline of the project

- Steve Plimpton, Sandia National Laboratories (retired)
- Axel Kohlmeyer, Temple University
- Aidan Thompson, Sandia National Laboratories
 - Richard Berger, Los Alamos National Laboratory
- Germain Clavier, Laboratoire CIMAP, France
- Joel Clemmer, Sandia National Laboratories
- Jacob R. Gissinger, Stevens Institute of Technology
- James Goff, Sandia National Laboratories
- Meg McCarthy, Sandia National Laboratories
- Stan Moore, Sandia National Laboratories
- Trung Nguyen, University of Chicago

These eleven core developers span six institutions, with no single
institution employing more than half.

### 16. Please describe the project's decision-making process

Technical and administrative decisions are made by achieving broad
consensus across the full group of core developers listed in Section 15.
This occurs primarily during monthly videoconference meetings, with
additional discussion and decision-making during annual in-person
developer meetings and on GitHub.  The three lead developers (Section
14) facilitate this process and act as tie-breakers of last resort, but
do not constitute a separate governing body and do not overrule core
developer consensus.  These practices are being formalized in the
governance document described in Section 7.

### 17. What is the maturity level of your project?

LAMMPS is over 30 years old and has existed in its current form — the
C++ code base distributed under the GPL — for more than 20 years.  Over
that time it has undergone sustained growth in functionality, platform
support, users, and stakeholders: the code base has grown to about 1.4
million lines of code in roughly 7,000 source files, organized into a
core plus 94 optional packages, with over 50,000 commits from more than
500 contributors, and a sustained rate of roughly 3,000–4,000 commits
per year over the past five years.

### 18. Please list the project's official communication channels

- Email: developers@lammps.org
- Slack: lammps.slack.com
- User forum: https://matsci.org/c/lammps

The Slack workspace is reserved for communication between developers and access is
on request only to steer all user support inquiries to the forum on
[MatSci.org](https://matsci.org/c/lammps).


### 19. Please list the project's social media accounts

YouTube: https://www.youtube.com/@templelammps

### 20. Please describe any existing financial sponsorships

Sandia National Laboratories supports Thompson, Moore, Clemmer, Goff,
and McCarthy via employment, and Kohlmeyer in part via a contract with
Temple University.  Temple University supports Kohlmeyer as well.
Individual core developers at other institutions are supported by their
employers and by research grants that include LAMMPS development
activities.

### 21. Please describe the project's infrastructure needs or requests

The LAMMPS project needs help managing both commmunity and infrastructure.
The most pressing need is sustainable hosting for its web
infrastructure.  The [lammps.org](https://www.lammps.org) website, the downloadable release
assets, and its static-analysis infrastructure are currently hosted at
Temple University and funded through research grants.  This arrangement
has the important property of not being encumbered by the security and
access restrictions of a national laboratory, but it is at risk of being
lost when the current research funding ends.

An additional need is testing capacity.  LAMMPS is a large
project with many configuration permutations and supported platforms
(see the testing dashboard at
[lammps-test-results](https://lammps.github.io/lammps-test-results/)), and the capacity limits
of GitHub's free tier for open source projects constrain how much of
this matrix can be tested, and how quickly.  Additional CI resources
would allow broader coverage, in particular for the growing regression
test suite.

Finally, LAMMPS needs support for managing meetings, including both developer
meetings and user training events.
