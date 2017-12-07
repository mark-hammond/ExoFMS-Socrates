# EXOFMS-SOC

![socrates](misc/socrates_structure.001.png)


This is an unfinished repo for the coupling of Socrates to Exo-FMS

## Installation

Download the latest Socrates version from XXX.

Install it to your srcmods folder.

Copy the XX files from this repo.

## Usage

 
```fortran
    id_soc_olr = &
    register_diag_field ( soc_mod_name, 'soc_olr', axes(1:2), Time, &
               'outgoing longwave radiation', &
               'watts/m2', missing_value=missing_value               )
```

Markdown | Less | Pretty
--- | --- | ---
*Still* | `renders` | **nicely**
1 | 2 | 3

Set the spectral file at XX, from:

1) XX
2) XX
3) XX

Set the active gases in XX.

Set a constant mixing ratio in XX.

Couple a dynamic mixing ratio in XX

## Examples

Run the test XX for an Earthlike planet.

Run the test XX for a lava planet.

## Generating Spectral Files

See XX for detailed instructions.

##Acknowledgements

Thanks to the Socrates team and UK Met Office for developing and sharing the model.

Thanks to James Manners for his particular help in getting Socrates coupled to ExoFMS, and for suggesting it in the first place!
