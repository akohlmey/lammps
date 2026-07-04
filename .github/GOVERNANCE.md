# LAMMPS Project Governance

This document describes how the LAMMPS project is led, how technical and
administrative decisions are made, and how contributors can become part of
that process.  It formalizes practices the project has followed informally
for many years.

## 1. Roles

### 1.1 Contributors

Anyone who submits an issue or pull request, or takes part in discussions on
the project's GitHub repository or official communication channels
(developers@lammps.org, lammps.slack.com, https://matsci.org/c/lammps) is a
Contributor.  No approval or affiliation is required.

### 1.2 Core Developers

Core Developers are Contributors who have been granted commit access to the
mainline repository.  They review and merge pull requests, participate in
technical and administrative decisions, and attend the project's regular
developer meetings.

**Becoming a Core Developer:** A Contributor is nominated by an existing
Core Developer on the basis of a sustained record of high-quality
contributions and sound judgment in review.  The nomination is confirmed by
consensus of the existing Core Developers (see Section 3).

**Stepping down / Emeritus status:** Core Developers may step down
voluntarily at any time.  A Core Developer who has been inactive for an
extended period (nominally 18 months) may be moved to Emeritus status by
consensus of the remaining Core Developers; commit access is removed, but
Emeritus developers remain credited for their contributions and may request
reinstatement.

**Removal for cause:** Commit access may be revoked by consensus of the
Lead Developers in response to a substantiated Code of Conduct violation,
following the enforcement process of the project's
[Code of Conduct](CODE_OF_CONDUCT.md).

### 1.3 Collaborators

Contributor submitting pull requests with significant new features are
invited to become Collaborators on the project. They are added to the
`.github/CODEOWNERS.md` file for their files / packages to be automatically
flagged as reviewer of any changes to their contributed files and are
given Triage permissions to review any changes to their contribution.

### 1.4 Lead Developers

The Lead Developers are Core Developers who additionally hold overall
technical and administrative responsibility for the project: they
facilitate consensus-building, act as tie-breakers of last resort when
consensus cannot be reached among the Core Developers, sign off on
releases, and represent LAMMPS externally.

Lead Developers are drawn from, and remain, Core Developers; they do not
constitute a separate governing body and do not have authority to
overrule a clear Core Developer consensus.  New Lead Developers are
selected by the current Lead Developers from among the Core Developers,
in consultation with the Core Developer group.

## 2. Institutional diversity

The Core Developer group is intentionally composed of individuals
affiliated with multiple independent institutions.  The project is
committed to maintaining this diversity: no single institution may
employ more than half of the Core Developers, and no single institution
may control project decisions.

## 3. Decision-making process

Technical and administrative decisions — including feature and package
proposals, API changes, deprecations, committer nominations, and
governance changes — are made by seeking broad consensus among the Core
Developers.  This happens primarily through:

- discussion and review on GitHub pull requests and issues;
- monthly developer videoconference meetings;
- an annual in-person developer meeting, used for larger architectural and
  strategic decisions.

If consensus cannot be reached through discussion, the Lead Developers
make the final decision, informed by the range of views expressed by the
Core Developers.

## 4. Release process

The Lead Developers decide when feature releases are cut and when a
release is promoted to a stable release, and are responsible for
publishing and announcing releases.  The release mechanics — branch
structure, tagging conventions, testing requirements, and integrity
protection — are documented in the LAMMPS manual and in
`.github/release_steps.md`.

## 5. Code of Conduct

All participants in the project are expected to follow the [LAMMPS Code
of Conduct](CODE_OF_CONDUCT.md).  Violations may be reported
confidentially to the project team at developers@lammps.org; reports are
reviewed and investigated by the Lead Developers as described in the
Code of Conduct's enforcement section.

## 6. Amending this document

Changes to this governance document are proposed as pull requests and
adopted by consensus of the Core Developers, following the same process
as other administrative decisions (Section 3).
