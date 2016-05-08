define(['app'], function (App) {

    App.Model.Ui = App.Model.Ui || {};

    App.Model.Ui.FilterEntry = App.Backbone.Model.extend({
        defaults: {
            active: false,
            name: null,
            count: 0
        }
    });

    App.Model.Ui.FilterEntries = App.Backbone.Collection.extend({
        model: App.Model.Ui.FilterEntry
    });

    App.Model.Ui.FilterGroup = App.Backbone.Model.extend({
        defaults: {
            active: false,
            name: null,
            label: null,
            entries: new App.Model.Ui.FilterEntries()
        }
    });

    App.Model.Ui.FilterGroups = App.Backbone.Collection.extend({
        model: App.Model.Ui.FilterGroup
    });

    return App.Model.Ui.FilterGroups;
});
