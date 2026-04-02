function mysave(fig, name)
    change_to_dark = strcmp(fig.Theme.BaseColorStyle, "dark");
    if change_to_dark
        theme(fig,"light");
    end
    saveas(fig, append(name, '.png'));
    saveas(fig, append(name, '.fig'));
    if change_to_dark
        theme(fig,"dark");
    end
end