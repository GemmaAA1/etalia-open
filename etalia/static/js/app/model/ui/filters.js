define(['app'], function (App) {

    App.Model.Ui = App.Model.Ui || {};

    App.Model.Ui.FilterEntry = App.Backbone.RelationalModel.extend({
        defaults: {
            active: false,
            label: null,
            count: 0
        }
    });

    App.Model.Ui.FilterEntries = App.Backbone.Collection.extend({
        model: App.Model.Ui.FilterEntry
    });

    App.Model.Ui.FilterGroup = App.Backbone.RelationalModel.extend({
        idAttribute: 'name',

        defaults: {
            active: false,
            name: null,
            label: null,
            entries: new App.Model.Ui.FilterEntries()
        },

        relations: [{
            type: Backbone.HasMany,
            key: 'entries',
            relatedModel: App.Model.Ui.FilterEntry,
            collectionType: App.Model.Ui.FilterEntries
        }]
    });

    App.Model.Ui.FilterGroups = App.Backbone.Collection.extend({
        model: App.Model.Ui.FilterGroup,

        parse: function(response) {
            return response.groups;
        }
    });

    return App.Model.Ui.FilterGroups;
});
