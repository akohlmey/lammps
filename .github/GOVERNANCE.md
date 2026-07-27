---
title: "Governance"
description: "How the LAMMPS project is managed"
---

<sub>(this page and its content is currently a work in progress and not final)</sub>

<!-- # Governance of the LAMMPS project -->

All information in this file is pursuant to the LAMMPS Project
[Technical Charter](/about/charter).  In case of any conflict or
difference in language, the LAMMPS Project Technical Charter takes
priority over this document.

## Technical Steering Committee

### Role

The role of the Technical Steering Committee (TSC) is to provide
technical direction to the project.  The TSC will vote on any matters on
which the community is unable to reach consensus.

### Membership

Members can be added to the TSC by a majority vote of the TSC.  Members
may be removed from the TSC by a 2/3 vote of the TSC.  If a TSC member
has been inactive for over 6 months, the TSC must hold a vote on whether
to remove that member from the TSC.

Current Membership:

1. Steve Plimpton  ([@sjplimp](https://github.com/sjplimp)), unaffiliated
1. Axel Kohlmeyer ([akohlmey](https://github.com/akohlmey)), Temple U.
1. Aidan Thompson ([@athomps](https://github.com/athomps)), SNL
1. Richard Berger ([@rbberger](https://github.com/rbberger)), LANL
1. Stan Moore ([@stanmoore1](https://github.com/stanmoore1)), SNL
1. Joel Clemmer, ([@jtclemm](https://github.com/jtclemm)), SNL
1. Jake Gissinger ([@jrgissing](https://github.com/jrgissing)), Stevens IT
1. Larry Fried (), LLNL
1. [TBD], LANL, ORNL,ANL, NVIDIA, Other stakeholders

### TSC Chair

The TSC will elect a chair.  The TSC chair runs TSC meetings and may
make interim decisions on urgent matters on behalf of the TSC, which may
be reviewed by the TSC at its next meeting.

Current chair: Aidan Thompson ([@athomps](https://github.com/athomps)), SNL

### Meetings and Notes

The TSC meets are held quarterly.  These meetings are open to the
public, and are held virtually.

Meeting notes can be found in [tsc-meeting-notes](tsc-meeting-notes)
in this repository.

## Other Public LAMMPS fora

### Monthly Developer Meetings

Public LAMMPS developer meetings are held on the second Tuesday of every
month that is not a TSC meeting, at 12 noon Albuquerque Time. These
meetings are held virtually over
[Teams](https://teams.microsoft.com/l/meetup-join/19%3ameeting_YWVlZjUwYzYtMTBhOS00MjM2LTgzMjgtOWU1NzdkNTUwZDZh%40thread.v2/0?context=%7b%22Tid%22%3a%227ccb5a20-a303-498c-b0c1-29007381b574%22%2c%22Oid%22%3a%22e1d7203f-43a9-4df0-a5fe-cbae01e46797%22%7d).

### LAMMPS Slack

The LAMMPS project uses slack for day-to-day developer communication at
[lammps.slack.com](https://lammps.slack.com).  Invitations to the LAMMPS
slack instance can be obtained by sending an email to slack@lammps.org.

<!-- https://lammps.slack.com/archives/C1G0VMHSL -->

### LAMMPS Github

LAMMPS [issues](https://github.com/lammps/lammps/issues) and [pull
requests](https://github.com/lammps/lammps/pulls), are tracked on the
[LAMMPS GitHub page](https://github.com/lammps/lammps).

### Announcements

Announcements will be posted on the LAMMPS website
[lammps.org](https://www.lammps.org/) and on
[MatSci.org/lammps](https://matsci.org/lammps).

## Teams within the LAMMPS project

### Administrative Access to the Github Repository

Only members of the TSC may have administrative access to the LAMMPS
repository. TSC members will be given access as needed. This list is not
generally published for operational security reasons.

### Other teams within the LAMMPS project

There are other teams within the LAMMPS GitHub project which control
access to various permissions on the repository. These teams are:

- `tsc`: see above
- `maintainers`: can merge an approved pull request to LAMMPS
- `core-developers`: can override a request for change, effectively
  allowing them to approve a pull request

See the relevant teams in the LAMMPS GitHub repository for membership
of these teams.  Membership in these teams is controlled by the TSC.

-----

A [copy of this document](https://github.com/akohlmey/lammps/blob/internal-docs/.github/GOVERNANCE.md)
is also available in the LAMMPS git repository as `.github/GOVERNANCE.md`.
