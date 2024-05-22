# Gas mixtures

Two pet peeves I'm constantly having from the codebase are:
- entropy
- lack of adiabatic expansion.

## Entropy

Yeah, I know that the comments honestly claim that

```
It's arguable whether this should even be called entropy anymore. It's more "based on" entropy than actually entropy now.
```

The git history helpfully shows the comment that is lost in the mists of time in the current code version -- that it is actually based on [Sackur-Tetrode equation](https://en.wikipedia.org/wiki/Sackur%E2%80%93Tetrode_equation) (for ideal monoatomic gas).

This warrants extra meditation on how the blasphemous unit of 
$\mathrm{\dfrac{mol^{7/3}\,K}{J^{2/3}\,kg^{2/3}\,\mathcal{l}}}$
appears in the current form.