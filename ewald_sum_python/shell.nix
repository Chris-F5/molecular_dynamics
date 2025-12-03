let
  nixpkgs = fetchTarball "https://github.com/NixOS/nixpkgs/tarball/nixos-24.05";
  pkgs = import nixpkgs { config = {}; overlays = []; };
in

pkgs.mkShellNoCC {
  packages = with pkgs; [
    pkgs.python312Full
    pkgs.python312Packages.numpy
    pkgs.python312Packages.scipy
    pkgs.python312Packages.matplotlib

    # tests
    pkgs.python312Packages.hypothesis
    pkgs.python312Packages.pytest

    # lsp
    pkgs.python312Packages.python-lsp-server
    pkgs.python312Packages.pylsp-rope
    pkgs.python312Packages.pyflakes
    pkgs.python312Packages.pycodestyle
  ];
}
