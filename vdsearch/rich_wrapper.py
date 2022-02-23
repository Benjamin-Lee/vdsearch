import rich_click as click
from rich_click import RichCommand as Command, RichGroup as Group
import typer

click.rich_click.USE_MARKDOWN = True
click.rich_click.MAX_WIDTH = 88

click.rich_click.COMMAND_GROUPS = {
    "vdsearch": [
        {"name": "Easy Workflows (recommended)", "commands": ["easy-search"]},
        {
            "name": "Download Datasets",
            "commands": ["download-cms", "download-viroiddb"],
        },
        {"name": "Internal (dangerous)", "commands": ["internal"]},
    ]
}
click.rich_click.STYLE_HELPTEXT = ""
click.rich_click.ERRORS_EPILOGUE = (
    "If you have any questions, please contact me at https://benjamindlee.com"
)
click.rich_click.FOOTER_TEXT = "BLWAG"


class MyTyper(typer.Typer):
    """A custom subclassed version of typer.Typer to allow rich help"""

    def __init__(self, *args, cls=Group, **kwargs) -> None:
        super().__init__(*args, cls=cls, **kwargs)

    def command(self, *args, cls=Command, **kwargs) -> typer.Typer.command:
        return super().command(*args, cls=cls, **kwargs)
