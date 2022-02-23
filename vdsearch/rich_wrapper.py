import rich_click as click
from rich_click import RichCommand as Command, RichGroup as Group
import typer
from click_didyoumean import DYMGroup

click.rich_click.USE_MARKDOWN = True
click.rich_click.MAX_WIDTH = 88

click.rich_click.COMMAND_GROUPS = {
    "vdsearch": [
        {"name": "Easy Workflows (recommended)", "commands": ["easy-search"]},
        {
            "name": "Download Datasets",
            "commands": ["download-cms", "download-viroiddb"],
        },
        {
            "name": "Search Steps (advanced)",
            "commands": [
                "find-circs",
                "dedup",
                "cluster",
            ],
        },
        {"name": "Internal (dangerous)", "commands": ["internal"]},
    ]
}
click.rich_click.STYLE_HELPTEXT = ""
click.rich_click.ERRORS_EPILOGUE = "[dim]If you have any questions or need help, please contact me at [cyan]https://benjamindlee.com"


class MyRichGroup(Group, DYMGroup):
    """Wrapper to enable Did You Mean"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class MyTyper(typer.Typer):
    """A custom subclassed version of typer.Typer to allow rich help"""

    def __init__(self, *args, cls=MyRichGroup, **kwargs) -> None:
        super().__init__(*args, cls=cls, **kwargs)

    def command(self, *args, cls=Command, **kwargs) -> typer.Typer.command:
        return super().command(*args, cls=cls, **kwargs)
