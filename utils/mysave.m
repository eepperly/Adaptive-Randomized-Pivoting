function mysave(fig, name)
    change_to_dark = strcmp(fig.Theme.BaseColorStyle, "dark");
    if change_to_dark
        theme(fig,"light");
    end
    saveas(fig, [name '.png']);
    saveas(fig, [name '.fig']);
    if change_to_dark
        theme(fig,"dark");
    end
end